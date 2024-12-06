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
//    LOG << "For node "<< hn << " in cluster " << from;
    // we want to compute the gain of moving the node out of the current cluster and
    // in to a neighboring cluster

    // map to store incident cluster IDs and the gain in pi modularity associated
    // with them
    std::unordered_map<PartitionID, double> delta_supt;
    double theta = 0.3;

    // iterate over all incident edges of hn to compute change in support of incident hyperedges
    // if the hypernode is moved to the corresponding cluster
    for (const HyperedgeID& he : phg.incidentEdges(hn)) {
        //LOG << "For hn " << hn << " and for he " << he <<":";
        // map to store cluster ID as key and fraction of pins in that cluster as loyalty value
        std::unordered_map<PartitionID, double> per_cluster_loyalty;
        double totalEdgeWeight = 0;

        // get total number of pins of the hyperedge
        //HypernodeID totalPins = phg.edgeSize(he);

        // go over all pins of the hyperedge and populate the map with loyalties for each incident cluster
        PartitionID clusterID = 0;
        for (const HypernodeID &pin: phg.pins(he)) {
            if(pin != hn) {
                clusterID = phg.partID(pin);
                //per_cluster_loyalty[clusterID] += static_cast<double>(phg.nodeWeight(pin));
                per_cluster_loyalty[clusterID] += (static_cast<double>(phg.nodeWeight(pin))/static_cast<double>(phg.edgeSize(he)));
            }
            //totalEdgeWeight += phg.nodeWeight(pin);
            //totalEdgeWeight += phg.nodeStrength(pin);
            totalEdgeWeight += (static_cast<double>(phg.nodeWeight(pin))/static_cast<double>(phg.edgeSize(he)));
        }

        // loyalty of hyperedge if hn is in its own cluster
        //double l_1 = static_cast<double>(phg.nodeWeight(hn)) / totalEdgeWeight;
        double l_1 = (static_cast<double>(phg.nodeWeight(hn))/static_cast<double>(phg.edgeSize(he))) / totalEdgeWeight;
//        LOG << "l_1 = " << (static_cast<double>(phg.nodeWeight(hn))/static_cast<double>(phg.edgeSize(he))) << " / " << totalEdgeWeight << " = " << l_1;
        double l_1_rho = compute_loyalty_rho(l_1, totalEdgeWeight, theta);

        for (const auto &pair: per_cluster_loyalty) {
            clusterID = pair.first;

            // loyalty of hyperedge if hn is kept in its current cluster
            double l_2 = pair.second / totalEdgeWeight;
//            LOG << "l_2 (" << clusterID << ") = " << pair.second << " / " << totalEdgeWeight << " = " << l_2;
            double l_2_rho = compute_loyalty_rho(l_2, totalEdgeWeight, theta);

            // Process the partition ID and loyalty value

            // loyalty of hyperedge if hn is sent to current clusterID
            double l_3 = l_1 + l_2;
            double l_3_rho = compute_loyalty_rho(l_3, totalEdgeWeight, theta);

            delta_supt[clusterID] += (l_3_rho - l_1_rho - l_2_rho);
//            LOG << "delta supt to = " << clusterID << " is " << delta_supt[clusterID];
        }
    }

    // compute pi modularity change if the hypernode is removed from its current cluster
    double change_in_pi_modularity_u_from_C = 0;

    // if hn is the only node in the cluster, removing it does not change pi mod score
    if(phg.partWeight(from) != phg.nodeWeight(hn)) {
        change_in_pi_modularity_u_from_C = deltaPIRemove(phg, hn, from, delta_supt[from]);
    }

    //constexpr Gain max_Gain = std::numeric_limits<Gain>::max();

    for (const auto& pair : delta_supt) {
      PartitionID clusterID = pair.first;
      if (clusterID != phg.partID(hn)) {
          double delta_supt_C = pair.second;
          double change_in_pi_modularity_u_to_C = deltaPI(phg, hn, clusterID, delta_supt_C);
          double net_change_in_pi_modularity = change_in_pi_modularity_u_to_C + change_in_pi_modularity_u_from_C;
//          LOG << "Node " << hn << " has strength = " << phg.nodeStrength(hn) << " and cluster " << clusterID
//          << " has strength = " << phg.partVolume(clusterID);
          tmp_scores[clusterID] = static_cast<HyperedgeWeight>(std::floor(net_change_in_pi_modularity * 10000000));
//          LOG << "Node "<< hn << " to cluster " << clusterID << " has pi_mod gain " << net_change_in_pi_modularity
//                      << " with from = " << change_in_pi_modularity_u_from_C <<
//                      " and to = " << change_in_pi_modularity_u_to_C;;
      }
    }
    //LOG << "---------";
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
                 const HypernodeID hn,
                 PartitionID new_cluster,
                 double delta_supt_C) {
      // this function returns the change in modularity on moving hn to new_cluster
//      LOG << "Computing deltaPI for node " << hn << " to cluster " << new_cluster;
      auto vol_H = static_cast<double>(phg.topLevelTotalVertexDegree());
      auto m = static_cast<double>(phg.topLevelNumEdges());

      // Calculate gamma
      const double gamma = (vol_H - 2 * m) / (vol_H - m);
      double theta = 0.3;

      // volume of new_cluster
      double vol_C = phg.partVolume(new_cluster);
      double eta_C = theta * (1.0 - (vol_C / m));

      // volume of cluster containing only hn
      double vol_hn = phg.nodeStrength(hn);
      double eta_hn = theta * (1.0 - (vol_hn / m));

      // volume of new_cluster with hn
      double vol_C_with_hn = vol_C + vol_hn;
      double eta_C_with_hn = theta * (1.0 - (vol_C_with_hn / m));

//      LOG << "Vol_hn = " << vol_hn << "; Vol_C = " << vol_C << "; Vol_C_with_hn = " << vol_C_with_hn;
//      LOG << "Delta_hn = " << expected_edges_in_cluster(gamma,eta_hn) << "; Delta_C = " << expected_edges_in_cluster(gamma,eta_C) << "; Delta_C_with_hn = " << expected_edges_in_cluster(gamma, eta_C_with_hn);
//      LOG << "delta_supt / m = " << delta_supt_C << " / " << m << " = " << delta_supt_C/m;

      double change_in_expected_edges = (delta_supt_C / m) + expected_edges_in_cluster(gamma, eta_C) +
              expected_edges_in_cluster(gamma, eta_hn) - expected_edges_in_cluster(gamma, eta_C_with_hn);

//      LOG << "deltaPiTo " << new_cluster << " = " << change_in_expected_edges;
      return change_in_expected_edges;
    }

    template<typename PartitionedHypergraph>
    double deltaPIRemove(const PartitionedHypergraph& phg,
                   const HypernodeID hn,
                   PartitionID old_cluster,
                   double delta_supt_C) {
        // this function returns the change in modularity on moving hn to new_cluster
//        LOG << "Computing deltaPI for node " << hn << " out from cluster " << old_cluster;
        auto vol_H = static_cast<double>(phg.topLevelTotalVertexDegree());
        auto m = static_cast<double>(phg.topLevelNumEdges());

        // Calculate gamma
        const double gamma = (vol_H - 2 * m) / (vol_H - m);
        double theta = 0.3;

        // volume of old_cluster
        double vol_C = phg.partVolume(old_cluster);
        double eta_C = theta * (1.0 - (vol_C / m));
//        LOG << "For node " << hn << " with volume " << phg.nodeStrength(hn) << " to cluster " << old_cluster << " vol_C = " << vol_C;
        // volume of cluster containing only hn
        double vol_hn = phg.nodeStrength(hn);
        double eta_hn = theta * (1.0 - (vol_hn / m));

        // volume of old_cluster without hn
        double vol_C_without_hn = vol_C - vol_hn;
        double eta_C_without_hn = theta * (1.0 - (vol_C_without_hn / m));
//        LOG << "Vol_hn = " << vol_hn << "; Vol_C = " << vol_C << "; Vol_C_without_hn = " << vol_C_without_hn;
//        LOG << "Delta_hn = " << expected_edges_in_cluster(gamma,eta_hn) << "; Delta_C = " << expected_edges_in_cluster(gamma,eta_C) << "; Delta_C_without_hn = " << expected_edges_in_cluster(gamma, eta_C_without_hn);
//        LOG << "delta_supt / m = " << delta_supt_C << " / " << m << " = " << delta_supt_C/m;

        double change_in_expected_edges = (delta_supt_C / m) + expected_edges_in_cluster(gamma, eta_C_without_hn) +
                                          expected_edges_in_cluster(gamma, eta_hn) - expected_edges_in_cluster(gamma, eta_C);

//        LOG << "deltaPiFrom = " << change_in_expected_edges;
        return -1 * change_in_expected_edges;
    }

    static double expected_edges_in_cluster(const double gamma, double eta) {
        // return expected edges in cluster according to Random Hypergraph Expansion Model
        double exp_value = (1.0-gamma) * (std::pow(1.0 - eta, 2)/(1.0 - gamma + (gamma * eta)));
        return exp_value;
    }

    static double compute_loyalty_rho(double loyalty, HypernodeWeight  totalEdgeWeight, double threshold) {
        if (loyalty >= threshold) {
            //return loyalty / std::log2((1.0 / loyalty) + 1.0);
            double log = std::log((1.0/loyalty)+1) / std::log(2);
            double linear_log = loyalty * (1.0/log);
            return totalEdgeWeight * linear_log;
        }
        return 0.0;
    }

};

}  // namespace mt_kahypar
