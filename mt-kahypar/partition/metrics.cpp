/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "mt-kahypar/partition/metrics.h"

#include <algorithm>
#include <cmath>
#include <tests/partition/refinement/flow_refiner_mock.h>

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/partition/refinement/gains/pimod/pimod_attributed_gains.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar::metrics {
double deltaPI(const PartitionedHypergraph &hypergraph, HypernodeID hn,
               PartitionID id, double c);

namespace {
template <typename PartitionedHypergraph, Objective objective>
struct ObjectiveFunction {};

template <typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::cut> {
  Gain operator()(const PartitionedHypergraph &phg,
                  const HyperedgeID &he) const {
    return phg.connectivity(he) > 1 ? phg.edgeWeight(he) : 0;
  }
};

template <typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::km1> {
  Gain operator()(const PartitionedHypergraph &phg,
                  const HyperedgeID &he) const {
    return std::max(phg.connectivity(he) - 1, 0) * phg.edgeWeight(he);
  }
};

template <typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::soed> {
  Gain operator()(const PartitionedHypergraph &phg,
                  const HyperedgeID &he) const {
    const PartitionID connectivity = phg.connectivity(he);
    return connectivity > 1 ? connectivity * phg.edgeWeight(he) : 0;
  }
};

template <typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::steiner_tree> {
  Gain operator()(const PartitionedHypergraph &phg,
                  const HyperedgeID &he) const {
    ASSERT(phg.hasTargetGraph());
    const TargetGraph *target_graph = phg.targetGraph();
    const Gain distance =
        target_graph->distance(phg.shallowCopyOfConnectivitySet(he));
    return distance * phg.edgeWeight(he);
  }
};

/* linear-over-log version used by the Pi-Mod paper */
inline double loyalty_rho(double l, double tw, double theta) {
  // if (l < theta)
  // return 0.0;
  // const double denom = std::log((1.0 / l) + 1.0) / std::log(2.0);
  // return l / denom; // --- FIX ---
  // tw = 1.0;
  if (l >= theta) {
    // return loyalty / std::log2((1.0 / loyalty) + 1.0);
    double log = std::log((1.0 / l) + 1) / std::log(2);
    double linear_log = l * (1.0 / log);
    return tw * linear_log;
    // return linear_log;
  }
  return 0.0;
}

template <typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::pimod> {
  Gain operator()(const PartitionedHypergraph &phg,
                  const HyperedgeID &he) const {
    // compute contribution of current hyperedge to hypergraph pi modularity
    double pi_mod_contribution = 0;

    // map to store cluster ID as key and fraction of pins in that cluster as
    // loyalty value std::unordered_map<PartitionID, double>
    // per_cluster_loyalty; vector of size = num. clusters adil: todo
    // potentially change for mem consumption
    std::vector<double> per_cluster_loyalty(phg.k(), 0);
    // LOG << "K = " << phg.k();

    // get total number of pins of the hyperedge
    // HypernodeID totalPins = phg.edgeSize(he);
    auto vol_H = static_cast<double>(phg.initialTotalVertexDegree());
    // auto vol_H = static_cast<double>(phg.topLevelTotalVertexDegree());
    auto m = static_cast<double>(phg.initialNumEdges());
    auto m_top = static_cast<double>(phg.topLevelNumEdges());
    // auto vol_H_top = static_cast<double>(phg.topLevelTotalVertexDegree());

    double theta = phg.getPiModTheta();
    double totalEdgeWeight = 0;

    // Calculate gamma
    double gamma = (vol_H - 2 * m) / (vol_H - m);

    // go over all pins of the hyperedge and populate the map with loyalties for
    // each incident cluster
    size_t pin_idx = 0;
    for (const HypernodeID &pin : phg.pins(he)) {
      PartitionID clusterID = phg.partID(pin);
      per_cluster_loyalty[clusterID] += phg.getNodeStrength(pin_idx, he);
      // per_cluster_loyalty[clusterID] +=
          // static_cast<double>(phg.nodeWeight(pin));
      totalEdgeWeight += phg.getNodeStrength(pin_idx, he);
      // totalEdgeWeight += static_cast<double>(phg.nodeWeight(pin));
      pin_idx++;
    }

    for (size_t cluster = 0; cluster < per_cluster_loyalty.size(); cluster++) {
      if (phg.partWeight(cluster) !=
          0) { // check if cluster is non-empty, i.e., it exists
        // Calculate eta for the current cluster C
        double vol_C = phg.partStrength(cluster);
        // auto vol_C = static_cast<double>(phg.partVolume(cluster));
        double eta = theta * (1.0 - (vol_C / m_top));
        // double eta = theta * (1.0 - (vol_C / vol_H_top));
        // Calculate expected edges in cluster according to Random Hypergraph
        // Expansion Model double expected_edges = (std::pow(1.0 - eta, 2) *
        // std::pow((1.0 + ((gamma * eta) / (1.0 - gamma))),-1))/m;
        double expected_edges = m_top* 
            std::pow((1.0 - eta), 2) *
             ((1.0 - gamma) / (1.0 - gamma + (gamma * eta)));

        if (per_cluster_loyalty[cluster] == 0) {
          // only remove exp edges to cluster / m from pi mod score
          pi_mod_contribution -= (expected_edges/m);
        } else if (per_cluster_loyalty[cluster] > 0) {
          const double loyalty_C =
              per_cluster_loyalty[cluster] / totalEdgeWeight;
          // Apply function rho to loyalty -> set to linear-over-log as
          // suggested by authors
          // double loyalty_rho = 0;
          // double log = 0;
          // double linear_log = 0;
          // if (loyalty_C >= theta) {
          //   // loyalty_rho = loyalty_C / std::log2((1.0 / loyalty_C) + 1.0);
          //   log = std::log((1.0 / loyalty_C) + 1) / std::log(2);
          //   linear_log = loyalty_C * (1.0 / log);
          //   // loyalty_rho = totalEdgeWeight * linear_log;
          //   loyalty_rho = linear_log;
          // }
          // const double rho = loyalty_rho(loyalty_C, totalEdgeWeight, theta) /
          // totalEdgeWeight;
          const double rho = loyalty_rho(loyalty_C, totalEdgeWeight, theta);

          // pi_mod_contribution += ((loyalty_rho - expected_edges) / m);
          pi_mod_contribution += rho - (expected_edges/m);
        }
      }
    }
    // auto final =
    //     static_cast<HyperedgeWeight>(std::floor(pi_mod_contribution *
    //     10000000));
    //  auto final = pi_mod_contribution;
    return -1 * pi_mod_contribution;
  }
};

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

template <typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::hmod> {
  Gain operator()(const PartitionedHypergraph &phg,
                  const HyperedgeID &he) const {
    // compute contribution of current hyperedge to hypergraph modularity
    const double tau = std::numeric_limits<double>::infinity();
    const double resolution = 1.0; // set to impact of degree tax
    // auto vol_H = static_cast<double>(phg.initialTotalVertexDegree());
    auto vol_H = static_cast<double>(phg.topLevelTotalVertexDegree());
    auto d = phg.edgeSize(he);
    double w_e = static_cast<double>(phg.edgeWeight(he));
    std::unordered_map<PartitionID, HypernodeID> pinsInCluster;

    for (const HypernodeID &pin : phg.pins(he)) {
      PartitionID clusterID = phg.partID(pin);
      ++pinsInCluster[clusterID];
    }

    // (A) Edge-contribution EC(e)

    // Step 1: Identify most populer cluster among pins
    HypernodeID c_max = 0;
    for (const auto &[cid, cnt] : pinsInCluster)
      c_max = std::max(c_max, cnt);

    double wdc_edge = 0.0;
    if (c_max > d / 2) {
      if (tau == std::numeric_limits<double>::infinity()) {
        wdc_edge = (c_max == d) ? 1.0 : 0.0;
      } else {
        wdc_edge = std::pow(static_cast<double>(c_max) / d, tau);
      }
    }

    const double EC = w_e * wdc_edge;

    // (B) Degree-tax DT(e)
    double DT = 0.0;

    for (HypernodeID c = d / 2 + 1; c <= d; ++c) {
      double wdc = 0.0;
      if (tau == std::numeric_limits<double>::infinity()) {
        wdc = (c == d) ? 1.0 : 0.0;
      } else {
        wdc = std::pow(static_cast<double>(c_max) / d, tau);
      }

      if (wdc == 0.0)
        continue;

      // TO DO: Fix num. edges of given size
      for (PartitionID cid = 0; cid < phg.k(); ++cid) {
        if (phg.partWeight(cid) != 0) {
          double vol_C = phg.partVolume(cid);
          double p = vol_C / vol_H;
          double pmf = binomPMF(c, d, p);
          // boost::math::binomial_distribution<double,double> binom(d, p);
          // double pmf = boost::math::pdf(binom, c);
          DT += w_e * wdc * pmf;
        }
      }
    }
    double contrib = EC - (resolution * DT);
    // std::cout << contrib << "\n";
    return -1 * contrib;
    // return -1 * static_cast<HyperedgeWeight>(contrib * 100);
  }
};

template <typename PartitionedHypergraph>
double print_hmod_contributions_no_scale(const PartitionedHypergraph &phg) {
  const double tau = std::numeric_limits<double>::infinity();
  const double resolution = 1.0;
  // auto vol_H = static_cast<double>(phg.initialTotalVertexDegree());
    auto vol_H = static_cast<double>(phg.topLevelTotalVertexDegree());
  double W = 0.0;

  double total_contribution = 0.0;
  for (const HyperedgeID he : phg.edges()) {
    const HypernodeID d = phg.edgeSize(he);
    const double w_e = static_cast<double>(phg.edgeWeight(he));
    W += phg.edgeWeight(he);

    /* --- count pins per block ------------------------------------------------
     */
    std::unordered_map<PartitionID, HypernodeID> pins_in_block;
    for (const HypernodeID pin : phg.pins(he)) {
      ++pins_in_block[phg.partID(pin)];
    }

    /* --- (A) Edge-contribution EC(e) ----------------------------------------
     */
    HypernodeID c_max = 0;
    for (auto &&[_, c] : pins_in_block)
      c_max = std::max(c_max, c);

    double wdc_edge = 0.0;
    if (c_max > d / 2) {
      wdc_edge = (tau == std::numeric_limits<double>::infinity())
                     ? (c_max == d ? 1.0 : 0.0)
                     : std::pow(static_cast<double>(c_max) / d, tau);
    }
    const double EC = w_e * wdc_edge;

    /* --- (B) Degree-tax DT(e) ----------------------------------------------
     */
    double DT = 0.0;
    for (HypernodeID c = d / 2 + 1; c <= d; ++c) {
      double wdc = (tau == std::numeric_limits<double>::infinity())
                       ? (c == d ? 1.0 : 0.0)
                       : std::pow(static_cast<double>(c) / d, tau);
      if (wdc == 0.0)
        continue;

      for (PartitionID pid = 0; pid < phg.k(); ++pid) {
        if (phg.partWeight(pid) == 0)
          continue;

        const double p = static_cast<double>(phg.partVolume(pid)) / vol_H;
        const double pmf = binomPMF(c, d, p); // exact binomial PMF
        DT += w_e * wdc * pmf;
      }
    }

    const double contrib = EC - resolution * DT;
    // std::cout << "hyperedge " << he << "  →  hmod-contrib = "
    //           << contrib << std::endl;

    total_contribution += contrib;
  }
  std::cout << "Numerator = " << total_contribution << std::endl;
  std::cout << "Modularity = " << total_contribution / W << std::endl;
  return total_contribution;
}

template <typename PartitionedHypergraph>
double print_pimod_contributions_no_scale(const PartitionedHypergraph &phg) {
  const double theta = phg.getPiModTheta(); // ρ-threshold
  auto vol_H = static_cast<double>(phg.initialTotalVertexDegree());
  // auto vol_H = static_cast<double>(phg.topLevelTotalVertexDegree());
  auto m = static_cast<double>(phg.initialNumEdges());
  auto m_top = static_cast<double>(phg.topLevelNumEdges());

  double self_count = 0;
  for(const HyperedgeID he: phg.edges()) {
    if(phg.edgeSize(he) == 1) self_count++;
  }
 
  vol_H -= self_count;
  m -= self_count;
  /* γ from the paper (Random Hyper-graph Expansion Model) */
  double gamma = (vol_H - 2.0 * m) / (vol_H - m);
  // LOG << "Gamma = " << gamma;

  auto expected_edges = [&](double eta) {
    // (1-γ)·(1-η)² /(1-γ+γη)           (no 1/m factor yet!)
    return (1.0 - gamma) *
           (std::pow(1.0 - eta, 2.0) / (1.0 - gamma + gamma * eta));
  };

  double numerator = 0.0; // Σ_e,C (ρ-E)
  // double rho_ch = 0;

  for (const HyperedgeID he : phg.edges()) {
    // if(phg.edgeSize(he) == 1) LOG << "Self loop";
    if (!phg.edgeIsEnabled(he))
      continue;

    /* --- gather pin weights per block ------------------------------------ */
    std::vector<double> pins_in_part(phg.k(), 0.0);
    double total_w = 0.0;

    std::size_t pin_idx = 0;
    for (const HypernodeID pin : phg.pins(he)) {
      const PartitionID pid = phg.partID(pin);
      // LOG << "Pin " << pin << " is in part " << pid;
      const double w = phg.getNodeStrength(pin_idx, he);
      // auto w = phg.nodeWeight(pin);
      pins_in_part[pid] += w;
      total_w += w;
      ++pin_idx;
    }

    /* --- per block contribution ------------------------------------------ */
    for (PartitionID pid = 0; pid < phg.k(); ++pid) {
      if (phg.partWeight(pid) == 0)
        continue; // empty block

      const double vol_C = phg.partStrength(pid); // vol(C)
      // auto vol_C = phg.partVolume(pid);
      const double eta = theta * (1.0 - vol_C / m_top);

      // we do not multiply exp_edges by m because we are summing over all edges
      // and for each cluster we include exp_edges/m for each cluster
      const double exp_edges = expected_edges(eta); // E[e∈C]

      if (pins_in_part[pid] == 0.0) {
        // LOG << "Hyperedge " << he << " contributes -" << exp_edges/m << " to
        // block " << pid;
        numerator -= exp_edges; // ρ = 0
        // rho_ch += exp_edges/m;
      } else {
        const double l = pins_in_part[pid] / total_w;      // loyalty
        const double rho = loyalty_rho(l, total_w, theta); // ρ(l)
        // rho_ch += exp_edges/m;
        numerator += rho - (exp_edges);
        // LOG << "Hyperedge " << he << " contributes -" <<  (rho - exp_edges/m)
        // << " to block " << pid;
      }
    }
    // LOG <<"-----";
  }
  // LOG << "Change in Exp = " << rho_ch;

  std::cout << "NumeratorE = " << numerator << '\n'
            << "ModularityE = " << numerator / m_top << std::endl;

  /* ------------------------------------------------------------------ *
   * 1) SUPPORT: accumulate per cluster                                 *
   * ------------------------------------------------------------------ */
  std::vector<double> supt(phg.k(), 0.0);
  // double num_self_loops = 0;

  for (const HyperedgeID he : phg.edges()) {
    if (!phg.edgeIsEnabled(he)) continue;
    // if(phg.edgeSize(he) == 1) {
      // num_self_loops+=1.0;
      // LOG << RED << "sl" << WHITE;
    // }

    /* pin weights per part */
    std::vector<double> pins_in_part2(phg.k(), 0.0);
    double total_we = 0.0;

    std::size_t pin_idx = 0;
    for (const HypernodeID pin : phg.pins(he)) {
      const PartitionID pid = phg.partID(pin);
      const double we        = phg.getNodeStrength(pin_idx, he);   // same weight you used before
      // auto we = phg.nodeWeight(pin);   // same weight you used before
      pins_in_part2[pid] += we;
      total_we           += we;
      ++pin_idx;
    }

    for (PartitionID pid = 0; pid < phg.k(); ++pid) {
      if (pins_in_part2[pid] == 0.0) continue;      // ρ = 0 → contributes nothing to support
      const double l   = pins_in_part2[pid] / total_we;          // loyalty
      const double rho = loyalty_rho(l, total_we, theta);       // ρ(l)
      supt[pid] += rho;
    }
  }
  // LOG << "Num self = " << num_self_loops;

  // vol_H -= num_self_loops;
  // m-= num_self_loops;
  // gamma = (vol_H - 2.0 * m) / (vol_H - m);
  // LOG << "Updated gamma = " << gamma;

  /* ------------------------------------------------------------------ *
   * 2) EXPECTATION & MODULARITY: iterate over clusters                 *
   * ------------------------------------------------------------------ */
  double numeratorC = 0.0;
  // double total_vol = 0.0;

  for (PartitionID pid = 0; pid < phg.k(); ++pid) {
    if (phg.partWeight(pid) == 0) continue;       // skip empty clusters

    const double vol_C = phg.partStrength(pid);               // vol(C)
    // total_vol += vol_C;
    // const double vol_C = phg.partStrength(pid);               // vol(C)
    const double eta   = theta * (1.0 - vol_C / m_top);

    const double exp_edges = expected_edges(eta)*m_top;         // divide by m once, here
    numeratorC += (supt[pid] - exp_edges);
    // LOG << "Cluster " << pid << " gives " << supt[pid] << " - " << exp_edges << " = " << (supt[pid] - exp_edges)/m_top;
  }

  // LOG << "Vol_H = " << vol_H << " and Total_vol = " << total_vol;
  // LOG << "Alt vol_h = " << phg.topLevelTotalVertexDegree();
  // LOG << "M = " << m;
  LOG << "NumeratorC = " << numeratorC;
  LOG << "ModularityC = " << numeratorC / m_top;

  return numerator; // raw Σ(ρ-E)
}

template <Objective objective, typename PartitionedHypergraph>
Gain compute_objective_parallel(const PartitionedHypergraph &phg) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  // if constexpr (objective == Objective::hmod) {
  //   /*  discard the return value – we only want the side effect */
  //   (void)print_hmod_contributions_no_scale(phg);
  // }
  // if constexpr (objective == Objective::pimod) {
  //   (void)print_pimod_contributions_no_scale(phg);
  // }
  tbb::enumerable_thread_specific<Gain> obj(0);
  phg.doParallelForAllEdges(
      [&](const HyperedgeID he) { obj.local() += func(phg, he); });
  return obj.combine(std::plus<>()) / (PartitionedHypergraph::is_graph ? 2 : 1);
}

template <Objective objective, typename PartitionedHypergraph>
Gain compute_objective_sequentially(const PartitionedHypergraph &phg) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  Gain obj = 0;
  for (const HyperedgeID &he : phg.edges()) {
    obj += func(phg, he);
  }
  return obj / (PartitionedHypergraph::is_graph ? 2 : 1);
}

template <Objective objective, typename PartitionedHypergraph>
Gain contribution(const PartitionedHypergraph &phg, const HyperedgeID he) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  return func(phg, he);
}
} // namespace

template <typename PartitionedHypergraph>
Gain quality(const PartitionedHypergraph &hg, const Context &context,
             const bool parallel) {
  return quality(hg, context.partition.objective, parallel);
}

template <typename PartitionedHypergraph>
Gain quality(const PartitionedHypergraph &hg, const Objective objective,
             const bool parallel) {
  switch (objective) {
  case Objective::cut:
    return parallel ? compute_objective_parallel<Objective::cut>(hg)
                    : compute_objective_sequentially<Objective::cut>(hg);
  case Objective::km1:
    return parallel ? compute_objective_parallel<Objective::km1>(hg)
                    : compute_objective_sequentially<Objective::km1>(hg);
  case Objective::soed:
    return parallel ? compute_objective_parallel<Objective::soed>(hg)
                    : compute_objective_sequentially<Objective::soed>(hg);
  case Objective::steiner_tree:
    return parallel
               ? compute_objective_parallel<Objective::steiner_tree>(hg)
               : compute_objective_sequentially<Objective::steiner_tree>(hg);
  case Objective::pimod:
    return parallel ? compute_objective_parallel<Objective::pimod>(hg)
                    : compute_objective_sequentially<Objective::pimod>(hg);
  case Objective::hmod:
    return parallel ? compute_objective_parallel<Objective::hmod>(hg)
                    : compute_objective_sequentially<Objective::hmod>(hg);
  default:
    throw InvalidParameterException("Unknown Objective");
  }
  return 0;
}

template <typename PartitionedHypergraph>
Gain contribution(const PartitionedHypergraph &hg, const HyperedgeID he,
                  const Objective objective) {
  switch (objective) {
  case Objective::cut:
    return contribution<Objective::soed>(hg, he);
  case Objective::km1:
    return contribution<Objective::km1>(hg, he);
  case Objective::soed:
    return contribution<Objective::soed>(hg, he);
  case Objective::steiner_tree:
    return contribution<Objective::steiner_tree>(hg, he);
  case Objective::pimod:
    return contribution<Objective::pimod>(hg, he);
  case Objective::hmod:
    return contribution<Objective::hmod>(hg, he);
  default:
    throw InvalidParameterException("Unknown Objective");
  }
  return 0;
}

template <typename PartitionedHypergraph>
bool isBalanced(const PartitionedHypergraph &phg, const Context &context) {
  size_t num_empty_parts = 0;
  for (PartitionID i = 0; i < context.partition.k; ++i) {
    if (phg.partWeight(i) > context.partition.max_part_weights[i]) {
      return false;
    }
    if (phg.partWeight(i) == 0) {
      num_empty_parts++;
    }
  }
  return context.partition.preset_type == PresetType::large_k ||
         num_empty_parts <= phg.numRemovedHypernodes();
}

template <typename PartitionedHypergraph>
double imbalance(const PartitionedHypergraph &hypergraph,
                 const Context &context) {
  ASSERT(context.partition.perfect_balance_part_weights.size() ==
         (size_t)context.partition.k);

  double max_balance =
      (hypergraph.partWeight(0) /
       static_cast<double>(context.partition.perfect_balance_part_weights[0]));

  for (PartitionID i = 1; i < context.partition.k; ++i) {
    const double balance_i =
        (hypergraph.partWeight(i) /
         static_cast<double>(
             context.partition.perfect_balance_part_weights[i]));
    max_balance = std::max(max_balance, balance_i);
  }

  return max_balance - 1.0;
}

template <typename PartitionedHypergraph>
double
approximationFactorForProcessMapping(const PartitionedHypergraph &hypergraph,
                                     const Context &context) {
  if (!PartitionedHypergraph::is_graph) {
    tbb::enumerable_thread_specific<Gain> approx_factor(0);
    hypergraph.doParallelForAllEdges([&](const HyperedgeID &he) {
      const size_t connectivity = hypergraph.connectivity(he);
      approx_factor.local() +=
          connectivity <= context.mapping.max_steiner_tree_size ? 1 : 2;
    });
    return static_cast<double>(approx_factor.combine(std::plus<>())) /
           hypergraph.initialNumEdges();
  } else {
    return 1.0;
  }
}

namespace {
#define OBJECTIVE_1(X)                                                         \
  Gain quality(const X &hg, const Context &context, const bool parallel)
#define OBJECTIVE_2(X)                                                         \
  Gain quality(const X &hg, const Objective objective, const bool parallel)
#define CONTRIBUTION(X)                                                        \
  Gain contribution(const X &hg, const HyperedgeID he,                         \
                    const Objective objective)
#define IS_BALANCED(X) bool isBalanced(const X &phg, const Context &context)
#define IMBALANCE(X)                                                           \
  double imbalance(const X &hypergraph, const Context &context)
#define APPROX_FACTOR(X)                                                       \
  double approximationFactorForProcessMapping(const X &hypergraph,             \
                                              const Context &context)
} // namespace

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(OBJECTIVE_1)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(OBJECTIVE_2)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(CONTRIBUTION)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(IS_BALANCED)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(IMBALANCE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(APPROX_FACTOR)
} // namespace mt_kahypar::metrics
