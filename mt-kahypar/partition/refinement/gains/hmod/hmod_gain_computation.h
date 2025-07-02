/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 *all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include <vector>

#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/refinement/gains/gain_computation_base.h"
#include "mt-kahypar/partition/refinement/gains/hmod/hmod_attributed_gains.h"

namespace mt_kahypar {
class HModGainComputation
    : public GainComputationBase<HModGainComputation, HModAttributedGains> {
  using Base = GainComputationBase<HModGainComputation, HModAttributedGains>;
  using RatingMap = typename Base::RatingMap;

  static constexpr bool enable_heavy_assert = false;

public:
  HModGainComputation(const Context &context,
                      bool disable_randomization = false)
      : Base(context, disable_randomization) {}

  // ! Precomputes the gain to all adjacent blocks.
  // ! Conceptually, we compute the gain of moving the node to an non-adjacent
  // block ! and the gain to all adjacent blocks assuming the node is in an
  // isolated block. ! The gain of that node to a block to can then be computed
  // by ! 'isolated_block_gain - tmp_scores[to]' (see gain(...))
  template <typename PartitionedHypergraph>
  void precomputeGains(const PartitionedHypergraph &phg, const HypernodeID hn,
                       RatingMap &tmp_scores, Gain &isolated_block_gain,
                       const bool) {
    ASSERT(tmp_scores.size() == 0, "Rating map not empty");
    // currently assigned cluster of the hypernode hn
    const PartitionID from = phg.partID(hn);
    // const double vol_hn = phg.nodeDegree(hn);
    const double vol_hn = phg.nodeVolume(hn);
    // const double W_edges = static_cast<double>(phg.initialNumEdges());
    //  LOG << "For node "<< hn << " in cluster " << from;
    //  we want to compute the gain of moving the node out of the current
    //  cluster and in to a neighboring cluster

    double resolution = 1.0;

    // --- cache of degree-tax values  DT(vol)  -----------------------
    // static std::unordered_map<std::uint64_t, double> dt_cache;
    // auto DT = [&](double volume)->double {
    //   std::uint64_t key = static_cast<std::uint64_t>(volume);
    //   auto it = dt_cache.find(key);
    //   if (it != dt_cache.end()) return it->second;

    //   const double vol_ratio = volume / phg.initialTotalVertexDegree();
    //   double dt = 0.0;
    //   //for (const auto& [d, cnt_d] : phg.edgeSizeCounters()) {      // map:
    //   size→count for (HyperedgeID d = 2; d <= phg.maxEdgeSize(); ++d) {
    //     HyperedgeID cnt_d = phg.edgeSizeCount(d);
    //     for (HyperedgeID c = d/2 + 1; c <= d; ++c) {                       //
    //     strict majority
    //       const double wdc = (c == d) ? 1.0 : 0.0;                 //  τ = ∞
    //       if (!wdc) continue;
    //       const double pmf = binomPMF(c, d, vol_ratio);
    //       dt += cnt_d * wdc * pmf;
    //     }
    //   }
    //   dt = dt / W;
    //   return dt_cache[key] = dt;
    // };
    // static std::unordered_map<std::uint64_t, double> dt_cache;
    auto DT = [&](double V) -> double {
      // const std::uint64_t key = static_cast<std::uint64_t>(V);
      // auto it = dt_cache.find(key);
      // if (it != dt_cache.end())
      //   return it->second;

      const double p = V / phg.topLevelTotalVertexDegree(); // volume ratio
      double sum = 0.0;
      for (HyperedgeID d = 2; d <= phg.maxEdgeSize(); ++d) {
        const HyperedgeID cnt_d = phg.edgeSizeCount(d);
        if (!cnt_d)
          continue;
        // τ = ∞  ⇒  w_dc = 1  only for  c = d   →   BinomPMF(d,d,p) = p^d
        sum += cnt_d * std::pow(p, static_cast<int>(d));
      }
      return sum;
      // return dt_cache[key] = sum / W_edges; // normalised
    };

    // (1) Degree-tax loss
    const double vol_from = static_cast<double>(phg.partVolume(from));
    const double dt_loss = DT(vol_from) - DT(vol_from - vol_hn);

    // (2) Edge-contribution loss
    double ec_loss = 0.0;                            // same for all dest
    std::unordered_map<PartitionID, double> ec_gain; // per destination

    double weight = 0;
    for (const HyperedgeID he : phg.incidentEdges(hn)) {
      weight += static_cast<double>(phg.edgeWeight(he));
      const int d = static_cast<int>(phg.edgeSize(he));
      if (d <= 1) {
        continue; // ignore single-pin edges
      }
      const double w_e = static_cast<double>(phg.edgeWeight(he));

      /* count pins *via stored counters* */
      const int cnt_from = static_cast<int>(phg.pinCountInPart(he, from));

      /* --- EC loss for <from> ------------------------------------------- */
      if (cnt_from == d) // edge was entirely in <from>
        ec_loss += w_e;  // will become cut after the move

      /* --- EC gain for each neighbour ----------------------------------- */
      for (const PartitionID pid : phg.connectivitySet(he)) {
        if (pid == from)
          continue; // not a destination

        const int cnt_to = static_cast<int>(phg.pinCountInPart(he, pid));

        if (cnt_to + 1 == d)   // after the move edge becomes internal
          ec_gain[pid] += w_e; // gains w_e
        if (cnt_to == d)       // edge is already internal → loses w_e
          ec_gain[pid] -= w_e;
      }
    }
    if(phg.nodeVolume(hn) != weight) { 
      LOG << RED << "Volume Stored = " << phg.nodeVolume(hn) << " and computed = " << weight << WHITE;
    }

    /* -------------------- (3)  assemble ΔQ for each dest block ------------ */
    for (auto &&[to, ec_g] : ec_gain) {
      const double vol_to = phg.partVolume(to);
      const double dt_gain = DT(vol_to + vol_hn) - DT(vol_to);

      // ec_loss /= W;
      // ec_g /= W;

      const double deltaQ = (ec_g - ec_loss)                    // ΔEC
                            - resolution * (dt_gain - dt_loss); // −γ·ΔDT

      // store negative value (KaHyPar minimises) as scaled integer
      //tmp_scores[to] = static_cast<HyperedgeWeight>(std::llround(deltaQ * 100));
      tmp_scores[to] = deltaQ;
      // std::cout << "Moving to cluster " << to << " has gain = " << deltaQ <<
      // std::endl;
    }

    isolated_block_gain = 0;
  }

  Gain gain(const Gain to_score, const Gain isolated_block_gain) {
    return isolated_block_gain - to_score;
  }

  void changeNumberOfBlocksImpl(const PartitionID) {
    // Do nothing
  }

  inline double binomPMF(const HypernodeID k, const HypernodeID n,
                         const double p) {
    if (p == 0.0)
      return (k == 0) ? 1.0 : 0.0;
    if (p == 1.0)
      return (k == n) ? 1.0 : 0.0;

    const double logComb =
        std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
    return std::exp(logComb + k * std::log(p) + (n - k) * std::log(1.0 - p));
  }

  /* --------------------------------------------------------------------
   * Strict-majority weighting  wdc(d,c;τ)
   *   – τ  = 0       majority   (all c > d/2   → weight 1)
   *   – τ  = ∞       strict     (only c = d     → weight 1)
   *   – else         (c/d)^τ    for c > d/2
   * ------------------------------------------------------------------*/
  inline double wdc(const int c, const int d, const double tau) {
    if (c <= d / 2)
      return 0.0; // minority
    if (std::isinf(tau))
      return (c == d) ? 1.0 : 0.0;
    if (tau == 0.0)
      return 1.0;                                     // majority
    return std::pow(static_cast<double>(c) / d, tau); // generic τ
  }
};
} // namespace mt_kahypar
