/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/datastructures/hypergraph_common.h"

namespace mt_kahypar {
/**
 * Utility function to compute loyalty (rho).
 * @param loyalty Loyalty value
 * @param threshold Threshold for applying the rho function
 * @return Computed rho value
 */
inline double compute_loyalty_rho(double loyalty, double totalEdgeWeight, double threshold) {
  if (loyalty >= threshold) {
    // return loyalty / std::log2((1.0 / loyalty) + 1.0);
    double log = std::log((1.0 / loyalty) + 1) / std::log(2);
    double linear_log = loyalty * (1.0 / log);
    return totalEdgeWeight * linear_log;
  }
  return 0.0;
}

/**
 * After moving a node, we perform a synchronized update of the pin count values
 * for each incident hyperedge of the node based on which we then compute an
 * attributed gain value.
 */
struct PiModAttributedGains {
  static double expected_edges_in_cluster(const double gamma, double eta) {
    // return expected edges in cluster according to Random Hypergraph Expansion Model
    double exp_value = (1.0 - gamma) * (std::pow(1.0 - eta, 2) / (1.0 - gamma + (gamma * eta)));
    return exp_value;
  }

  static double deltaPIRemove(double delta_supt_C, const SynchronizedEdgeUpdate& sync_update) {
    // this function returns the change in modularity on moving hn to new_cluster

    auto vol_H = static_cast<double>(sync_update.vol_H);
    auto m = static_cast<double>(sync_update.m);

    // Calculate gamma
    const double gamma = (vol_H - 2 * m) / (vol_H - m);
    double theta = sync_update.theta;

    // volume of old_cluster (if hn in it)
    double vol_C = sync_update.vol_From + sync_update.hn_strength;
    double eta_C = theta * (1.0 - (vol_C / m));
    // LOG << "ATT: For node " << sync_update.hn << " with volume " << sync_update.hn_strength << " to cluster " << sync_update.from << " vol_C = " << vol_C;

    // volume of cluster containing only hn
    double vol_hn = sync_update.hn_strength;
    double eta_hn = theta * (1.0 - (vol_hn / m));

    // volume of old_cluster without hn
    double vol_C_without_hn = sync_update.vol_From;
    double eta_C_without_hn = theta * (1.0 - (vol_C_without_hn / m));

    double change_in_expected_edges = (delta_supt_C / m) + ((expected_edges_in_cluster(gamma, eta_C_without_hn) +
                                                             expected_edges_in_cluster(gamma, eta_hn) - expected_edges_in_cluster(gamma, eta_C)) / sync_update.hn_degree);

    return -1 * change_in_expected_edges;
  }

  static double deltaPI(double delta_supt_C, const SynchronizedEdgeUpdate& sync_update) {
    // this function returns the change in modularity on moving hn to new_cluster

    auto vol_H = static_cast<double>(sync_update.vol_H);
    auto m = static_cast<double>(sync_update.m);

    // Calculate gamma
    const double gamma = (vol_H - 2 * m) / (vol_H - m);
    double theta = sync_update.theta;

    // volume of new_cluster
    double vol_C = sync_update.vol_To - sync_update.hn_strength;
    double eta_C = theta * (1.0 - (vol_C / m));
    // LOG << "ATT: For node " << sync_update.hn << " with volume " << sync_update.hn_strength << " to cluster " << sync_update.to << " vol_C = " << vol_C;
    // LOG << "vol_C = " << vol_C << " and eta_C = " << eta_C;

    // volume of cluster containing only hn
    double vol_hn = sync_update.hn_strength;
    double eta_hn = theta * (1.0 - (vol_hn / m));

    // LOG << "vol_hn = " << vol_hn << " and eta_hn = " << eta_hn;

    // volume of new_cluster with hn
    double vol_C_with_hn = sync_update.vol_To;
    double eta_C_with_hn = theta * (1.0 - (vol_C_with_hn / m));

    // LOG << "vol_C_with_hn = " << vol_C_with_hn << " and eta_C_with_hn = " << eta_C_with_hn;

    double change_in_expected_edges = (delta_supt_C / m) + ((expected_edges_in_cluster(gamma, eta_C) +
                                                             expected_edges_in_cluster(gamma, eta_hn) - expected_edges_in_cluster(gamma, eta_C_with_hn)) / sync_update.hn_degree);

    return change_in_expected_edges;
  }

  static HyperedgeWeight gain(const SynchronizedEdgeUpdate& sync_update) {
    // here, we receive volumes of To and From, and loyalty to To and From with hn already in To
    double theta = sync_update.theta;
    // LOG << "ATT: For node " << sync_update.hn << " and edge " << sync_update.he;
    // LOG << "For he = " << sync_update.he;

    // loyalty of hyperedge if hn is in its own cluster
    double l_1 = sync_update.hn_loyalty / sync_update.edge_weight_from_nodes;
    // LOG << "ATT: l_1 = " << sync_update.hn_weight << " / " << sync_update.edge_weight_from_nodes << " = " << l_1;
    // double l_1 = (static_cast<double>(sync_update.hn_weight) / sync_update.edge_weight_from_nodes;
    double l_1_rho = compute_loyalty_rho(l_1, sync_update.edge_weight_from_nodes, theta);

    // compute loyalty of he to From cluster (does not include hn_weight)
    // double l_2_From = (static_cast<double>(sync_update.loyalty_towards_from_part)/static_cast<double>(sync_update.edge_size)) / static_cast<double>(sync_update.edge_weight_from_nodes);
    double l_2_From = sync_update.loyalty_towards_from_part / sync_update.edge_weight_from_nodes;
    double l_2_From_rho = compute_loyalty_rho(l_2_From, sync_update.edge_weight_from_nodes, theta);
    // LOG << "ATT: l_2_From = " << sync_update.loyalty_towards_from_part << " / " << sync_update.edge_weight_from_nodes << " = " << l_2_From;

    // compute loyalty if hn is not in To yet
    // double l_2_To = (static_cast<double>(sync_update.loyalty_towards_to_part)/static_cast<double>(sync_update.edge_size)) / static_cast<double>(sync_update.edge_weight_from_nodes);
    double l_2_To = sync_update.loyalty_towards_to_part / sync_update.edge_weight_from_nodes;
    double l_2_To_rho = compute_loyalty_rho(l_2_To, sync_update.edge_weight_from_nodes, theta);
    // LOG << "ATT: l_2_To = " << sync_update.loyalty_towards_to_part << " / " << sync_update.edge_weight_from_nodes << " = " << l_2_To;

    // compute loyalty if hn is in From
    double l_3_From = l_2_From + l_1;
    double l_3_From_rho = compute_loyalty_rho(l_3_From, sync_update.edge_weight_from_nodes, theta);

    // compute loyalty if hn is in To
    double l_3_To = l_2_To + l_1;
    double l_3_To_rho = compute_loyalty_rho(l_3_To, sync_update.edge_weight_from_nodes, theta);

    double delta_supt_to_From_of_he = (l_3_From_rho - l_1_rho - l_2_From_rho);
    double delta_supt_to_C_of_he = (l_3_To_rho - l_1_rho - l_2_To_rho);
    // LOG << "ATT: delta supt from = " << sync_update.from << " is " << delta_supt_to_From_of_he;
    // LOG << "ATT: delta supt to = " << sync_update.to << " is " << delta_supt_to_C_of_he;

    // compute pi modularity change if the hypernode is removed from its cluster
    double change_in_pi_modularity_u_from_From = 0;
    if (sync_update.hn_weight != (sync_update.weight_From + sync_update.hn_weight)) {
      //   LOG << "HN: " << sync_update.hn_strength << " and vol From = " << sync_update.vol_From;
      //   LOG << RED << "ATT: Came" << WHITE;
      change_in_pi_modularity_u_from_From = deltaPIRemove(delta_supt_to_From_of_he, sync_update);
    }
    double change_in_pi_modularity_u_to_To = deltaPI(delta_supt_to_C_of_he, sync_update);

    double change_in_pi_modularity = change_in_pi_modularity_u_to_To + change_in_pi_modularity_u_from_From;
    // LOG << "ATT: Node " << sync_update.hn << " to cluster " << sync_update.to << " has pi_mod gain " <<
    //   change_in_pi_modularity << " with from = " << change_in_pi_modularity_u_from_From <<
    //   " and to = " << change_in_pi_modularity_u_to_To;
    auto final = static_cast<HyperedgeWeight>(std::floor(change_in_pi_modularity * 1000000));
    //auto final = change_in_pi_modularity;
    return -1* final;
  }
};
}  // namespace mt_kahypar
