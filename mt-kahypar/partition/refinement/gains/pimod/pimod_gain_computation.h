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
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
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

#include "mt-kahypar/partition/refinement/gains/gain_computation_base.h"
#include "mt-kahypar/partition/refinement/gains/pimod/pimod_attributed_gains.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

class PiModGainComputation : public GainComputationBase<PiModGainComputation, PiModAttributedGains> {
  using Base = GainComputationBase<PiModGainComputation, PiModAttributedGains>;
  using RatingMap = typename Base::RatingMap;

  static constexpr bool enable_heavy_assert = false;

 public:
    PiModGainComputation(const Context& context,
                     bool disable_randomization = false) :
    Base(context, disable_randomization) { }

  // ! Precomputes the gain to all adjacent blocks.
  // ! Conceptually, we compute the gain of moving the node to an non-adjacent block
  // ! and the gain to all adjacent blocks assuming the node is in an isolated block.
  // ! The gain of that node to a block to can then be computed by
  // ! 'isolated_block_gain - tmp_scores[to]' (see gain(...))
  template<typename PartitionedHypergraph>
  void precomputeGains(const PartitionedHypergraph& phg,
                       const HypernodeID hn,
                       RatingMap& tmp_scores,
                       Gain& isolated_block_gain,
                       const bool) {
    ASSERT(tmp_scores.size() == 0, "Rating map not empty");
    // currently assigned cluster of the hypernode hn
    PartitionID from = phg.partID(hn);

    // we want to compute the gain of moving the node out of the current cluster and
    // in to a neighboring cluster

    // map to store incident cluster IDs and the gain in pi modularity associated
    // with them
    std::unordered_map<PartitionID, double> delta_supt;

    // iterate over all incident edges of hn to compute change in support of incident hyperedges
    // if the hypernode is moved to the corresponding cluster
    for (const HyperedgeID& he : phg.incidentEdges(hn)) {
        LOG << "For he " << he <<":";
        // map to store cluster ID as key and fraction of pins in that cluster as loyalty value
        std::unordered_map<PartitionID, double> per_cluster_loyalty;

        // get total number of pins of the hyperedge
        HypernodeID totalPins = phg.edgeSize(he);

        // go over all pins of the hyperedge and populate the map with loyalties for each incident cluster
        PartitionID clusterID = 0;
        for (const HypernodeID &pin: phg.pins(he)) {
            if(pin != hn) {
                clusterID = phg.partID(pin);
                LOG << "Pin " << pin << " in cluster " << clusterID;
                per_cluster_loyalty[clusterID] += 1.0 / totalPins;
            }
        }

        // loyalty of hyperedge if hn is in its own cluster
        double l_1 = 1.0 / totalPins;
        double l_1_rho = l_1 / std::log2((1.0 / l_1) + 1.0);

        for (const auto &pair: per_cluster_loyalty) {
            clusterID = pair.first;

            // loyalty of hyperedge if hn is kept in its current cluster
            double l_2 = pair.second;
            double l_2_rho = l_2 / std::log2((1.0 / l_2) + 1.0);

            // Process the partition ID and loyalty value
            //std::cout << "Partition ID: " << clusterID << ", Loyalty: " << l_2 << std::endl;

            // loyalty of hyperedge if hn is sent to current clusterID
            double l_3 = l_1 + l_2;
            double l_3_rho = l_3 / std::log2((1.0 / l_3) + 1.0);

            delta_supt[clusterID] += (l_3_rho - l_1_rho - l_2_rho);
            LOG << "Delta supt for " << clusterID << " is now " << delta_supt[clusterID];
        }
    }

    // compute pi modularity change if the hypernode is removed from its current cluster
    double change_in_pi_modularity_u_from_C = deltaPIRemove(phg, from, delta_supt[from]);

    //constexpr Gain max_Gain = std::numeric_limits<Gain>::max();

    for (const auto& pair : delta_supt) {
      PartitionID clusterID = pair.first;
      double delta_supt_C = pair.second;
      double change_in_pi_modularity_u_to_C = deltaPI(phg, clusterID, delta_supt_C);
      double net_change_in_pi_modularity = change_in_pi_modularity_u_to_C + change_in_pi_modularity_u_from_C;
      //std::cout << deltaPI(phg, clusterID, delta_supt_C) << std::endl;
      tmp_scores[clusterID] =  (net_change_in_pi_modularity * 100000000000);
      LOG << "Cluster " << clusterID << " has pi_mod gain " << net_change_in_pi_modularity;
      //if (net_change_in_pi_modularity < -1.0 || net_change_in_pi_modularity > 1.0) {
      //  throw std::out_of_range("Value must be in the range [-1, 1]");
      //}
      // shift net_change_in_pi_modularity by 1.0 to make all values in the range [0,2]
      // map doubles in [0,2] to [0, 2^{32}-1] for unsigned int 32-bit ints
      //tmp_scores[clusterID] = static_cast<Gain>(std::round((net_change_in_pi_modularity+1.0) * (max_Gain/2.0)));
    }
    isolated_block_gain = 0;

  }

  HyperedgeWeight gain(const Gain to_score,
                       const Gain isolated_block_gain) {
    return isolated_block_gain - to_score;
  }

  void changeNumberOfBlocksImpl(const PartitionID) {
    // Do nothing
  }

  template<typename PartitionedHypergraph>
  double deltaPI(const PartitionedHypergraph& phg,
                 PartitionID new_cluster,
                 double delta_supt_C) {
      // this function returns the change in modularity on moving hn to new_cluster

      auto vol_H = static_cast<double>(phg.initialNumPins());
      auto m = static_cast<double>(phg.initialNumEdges());

      // Calculate gamma
      const double gamma = (vol_H - 2 * m) / (vol_H - m);
      double theta = 0.7;

      // volume of new_cluster
      auto vol_C = static_cast<double>(phg.partWeight(new_cluster));
      double eta_C = theta * (1.0 - (vol_C / vol_H));

      // volume of cluster containing only hn
      double vol_hn = 1.0;
      double eta_hn = theta * (1.0 - (vol_hn / vol_H));

      // volume of new_cluster with hn
      double vol_C_with_hn = vol_C + vol_hn;
      double eta_C_with_hn = theta * (1.0 - (vol_C_with_hn / vol_H));

      double change_in_expected_edges = (delta_supt_C / m) + expected_edges_in_cluster(gamma, eta_C) +
              expected_edges_in_cluster(gamma, eta_hn) - expected_edges_in_cluster(gamma, eta_C_with_hn);

      return change_in_expected_edges;
    }

    template<typename PartitionedHypergraph>
    double deltaPIRemove(const PartitionedHypergraph& phg,
                   PartitionID old_cluster,
                   double delta_supt_C) {
        // this function returns the change in modularity on moving hn to new_cluster

        auto vol_H = static_cast<double>(phg.initialNumPins());
        auto m = static_cast<double>(phg.initialNumEdges());

        // Calculate gamma
        const double gamma = (vol_H - 2 * m) / (vol_H - m);
        double theta = 0.7;

        // volume of old_cluster
        auto vol_C = static_cast<double>(phg.partWeight(old_cluster));
        double eta_C = theta * (1.0 - (vol_C / vol_H));

        // volume of cluster containing only hn
        double vol_hn = 1.0;
        double eta_hn = theta * (1.0 - (vol_hn / vol_H));

        // volume of old_cluster without hn
        double vol_C_without_hn = vol_C - vol_hn;
        double eta_C_without_hn = theta * (1.0 - (vol_C_without_hn / vol_H));

        double change_in_expected_edges = (delta_supt_C / m) + expected_edges_in_cluster(gamma, eta_C_without_hn) +
                                          expected_edges_in_cluster(gamma, eta_hn) - expected_edges_in_cluster(gamma, eta_C);

        return -1 * change_in_expected_edges;
    }

    static double expected_edges_in_cluster(const double gamma, double eta) {
        // return expected edges in cluster according to Random Hypergraph Expansion Model
        double exp_value = (1-gamma) * (std::pow(1.0 - eta, 2)/(1.0 - gamma + gamma * eta));
        return exp_value;
    }

};

}  // namespace mt_kahypar
