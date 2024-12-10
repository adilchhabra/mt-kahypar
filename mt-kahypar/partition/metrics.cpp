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

#include "mt-kahypar/partition/metrics.h"

#include <algorithm>
#include <cmath>
#include <tests/partition/refinement/flow_refiner_mock.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar::metrics {
double deltaPI(const PartitionedHypergraph& hypergraph, HypernodeID hn, PartitionID id, double c);

namespace {
template <typename PartitionedHypergraph, Objective objective>
struct ObjectiveFunction { };

template <typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::cut> {
  HyperedgeWeight operator() (const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    return phg.connectivity(he) > 1 ? phg.edgeWeight(he) : 0;
  }
};

template <typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::km1> {
  HyperedgeWeight operator() (const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    return std::max(phg.connectivity(he) - 1, 0) * phg.edgeWeight(he);
  }
};

template <typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::soed> {
  HyperedgeWeight operator() (const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    const PartitionID connectivity = phg.connectivity(he);
    return connectivity > 1 ? connectivity * phg.edgeWeight(he) : 0;
  }
};

template <typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::steiner_tree> {
  HyperedgeWeight operator() (const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    ASSERT(phg.hasTargetGraph());
    const TargetGraph* target_graph = phg.targetGraph();
    const HyperedgeWeight distance = target_graph->distance(phg.shallowCopyOfConnectivitySet(he));
    return distance * phg.edgeWeight(he);
  }
};

template <typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::pimod> {
  HyperedgeWeight operator() (const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    // compute contribution of current hyperedge to hypergraph pi modularity
    double pi_mod_contribution = 0;

    // map to store cluster ID as key and fraction of pins in that cluster as loyalty value
    // std::unordered_map<PartitionID, double> per_cluster_loyalty;
    // vector of size = num. clusters
    // adil: todo potentially change for mem consumption
    std::vector<double> per_cluster_loyalty(phg.k(), 0);
    // LOG << "K = " << phg.k();

    // get total number of pins of the hyperedge
    // HypernodeID totalPins = phg.edgeSize(he);
    auto vol_H = static_cast<double>(phg.initialTotalVertexDegree());
    auto m = static_cast<double>(phg.initialNumEdges());

    double theta = 0.5;
    double totalEdgeWeight = 0;

    // Calculate gamma
    double gamma = (vol_H - 2 * m) / (vol_H - m);

    // go over all pins of the hyperedge and populate the map with loyalties for each incident cluster
    size_t pin_idx = 0;
    for (const HypernodeID& pin : phg.pins(he)) {
      PartitionID clusterID = phg.partID(pin);
      per_cluster_loyalty[clusterID] += phg.getNodeStrength(pin_idx, he);
      totalEdgeWeight += phg.getNodeStrength(pin_idx, he);
      pin_idx++;
    }

    for (size_t cluster = 0; cluster < per_cluster_loyalty.size(); cluster++) {
      if (phg.partWeight(cluster) != 0) {       // check if cluster is non-empty, i.e., it exists
        // Calculate eta for the current cluster C
        double vol_C = phg.partVolume(cluster);
        double eta = theta * (1.0 - (vol_C / m));
        // Calculate expected edges in cluster according to Random Hypergraph Expansion Model
        // double expected_edges = (std::pow(1.0 - eta, 2) * std::pow((1.0 + ((gamma * eta) / (1.0 - gamma))),-1))/m;
        double expected_edges = (std::pow((1.0 - eta), 2) * ((1.0 - gamma) / (1.0 - gamma + (gamma * eta))));

        if (per_cluster_loyalty[cluster] == 0) {
          // only remove exp edges to cluster / m from pi mod score
          pi_mod_contribution -= (expected_edges / m);
        } else if (per_cluster_loyalty[cluster] > 0) {
          double loyalty_C = per_cluster_loyalty[cluster] / totalEdgeWeight;
          // Apply function rho to loyalty -> set to linear-over-log as suggested by authors
          double loyalty_rho = 0;
          double log = 0;
          double linear_log = 0;
          if (loyalty_C >= theta) {
            // loyalty_rho = loyalty_C / std::log2((1.0 / loyalty_C) + 1.0);
            log = std::log((1.0 / loyalty_C) + 1) / std::log(2);
            linear_log = loyalty_C * (1.0 / log);
            loyalty_rho = totalEdgeWeight * linear_log;
          }

          pi_mod_contribution += ((loyalty_rho - expected_edges) / m);
        }
      }
    }
    auto final = static_cast<HyperedgeWeight>(std::floor(pi_mod_contribution * 10000000));

    return -1* final;
  }
};

template <Objective objective, typename PartitionedHypergraph>
HyperedgeWeight compute_objective_parallel(const PartitionedHypergraph& phg) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  tbb::enumerable_thread_specific<HyperedgeWeight> obj(0);
  phg.doParallelForAllEdges([&](const HyperedgeID he) {
        obj.local() += func(phg, he);
      });
  return obj.combine(std::plus<>()) / (PartitionedHypergraph::is_graph ? 2 : 1);
}

template <Objective objective, typename PartitionedHypergraph>
HyperedgeWeight compute_objective_sequentially(const PartitionedHypergraph& phg) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  HyperedgeWeight obj = 0;
  for (const HyperedgeID& he : phg.edges()) {
    obj += func(phg, he);
  }
  return obj / (PartitionedHypergraph::is_graph ? 2 : 1);
}

template <Objective objective, typename PartitionedHypergraph>
HyperedgeWeight contribution(const PartitionedHypergraph& phg, const HyperedgeID he) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  return func(phg, he);
}
}

template <typename PartitionedHypergraph>
HyperedgeWeight quality(const PartitionedHypergraph& hg,
                        const Context& context,
                        const bool parallel) {
  return quality(hg, context.partition.objective, parallel);
}

template <typename PartitionedHypergraph>
HyperedgeWeight quality(const PartitionedHypergraph& hg,
                        const Objective objective,
                        const bool parallel) {
  switch (objective) {
    case Objective::cut:
      return parallel ? compute_objective_parallel<Objective::cut>(hg) :
             compute_objective_sequentially<Objective::cut>(hg);
    case Objective::km1:
      return parallel ? compute_objective_parallel<Objective::km1>(hg) :
             compute_objective_sequentially<Objective::km1>(hg);
    case Objective::soed:
      return parallel ? compute_objective_parallel<Objective::soed>(hg) :
             compute_objective_sequentially<Objective::soed>(hg);
    case Objective::steiner_tree:
      return parallel ? compute_objective_parallel<Objective::steiner_tree>(hg) :
             compute_objective_sequentially<Objective::steiner_tree>(hg);
    case Objective::pimod:
      return parallel ? compute_objective_parallel<Objective::pimod>(hg) :
             compute_objective_sequentially<Objective::pimod>(hg);
    default: throw InvalidParameterException("Unknown Objective");
  }
  return 0;
}

template <typename PartitionedHypergraph>
HyperedgeWeight contribution(const PartitionedHypergraph& hg,
                             const HyperedgeID he,
                             const Objective objective) {
  switch (objective) {
    case Objective::cut: return contribution<Objective::soed>(hg, he);
    case Objective::km1: return contribution<Objective::km1>(hg, he);
    case Objective::soed: return contribution<Objective::soed>(hg, he);
    case Objective::steiner_tree: return contribution<Objective::steiner_tree>(hg, he);
    case Objective::pimod: return contribution<Objective::pimod>(hg, he);
    default: throw InvalidParameterException("Unknown Objective");
  }
  return 0;
}

template <typename PartitionedHypergraph>
bool isBalanced(const PartitionedHypergraph& phg, const Context& context) {
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
double imbalance(const PartitionedHypergraph& hypergraph, const Context& context) {
  ASSERT(context.partition.perfect_balance_part_weights.size() == (size_t)context.partition.k);

  double max_balance = (hypergraph.partWeight(0) /
                        static_cast<double>(context.partition.perfect_balance_part_weights[0]));

  for (PartitionID i = 1; i < context.partition.k; ++i) {
    const double balance_i =
      (hypergraph.partWeight(i) /
       static_cast<double>(context.partition.perfect_balance_part_weights[i]));
    max_balance = std::max(max_balance, balance_i);
  }

  return max_balance - 1.0;
}

template <typename PartitionedHypergraph>
double approximationFactorForProcessMapping(const PartitionedHypergraph& hypergraph, const Context& context) {
  if (!PartitionedHypergraph::is_graph) {
    tbb::enumerable_thread_specific<HyperedgeWeight> approx_factor(0);
    hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
        const size_t connectivity = hypergraph.connectivity(he);
        approx_factor.local() += connectivity <= context.mapping.max_steiner_tree_size ? 1 : 2;
      });
    return static_cast<double>(approx_factor.combine(std::plus<>())) / hypergraph.initialNumEdges();
  } else {
    return 1.0;
  }
}

namespace {
#define OBJECTIVE_1(X) HyperedgeWeight quality(const X& hg, const Context& context, const bool parallel)
#define OBJECTIVE_2(X) HyperedgeWeight quality(const X& hg, const Objective objective, const bool parallel)
#define CONTRIBUTION(X) HyperedgeWeight contribution(const X& hg, const HyperedgeID he, const Objective objective)
#define IS_BALANCED(X) bool isBalanced(const X& phg, const Context& context)
#define IMBALANCE(X) double imbalance(const X& hypergraph, const Context& context)
#define APPROX_FACTOR(X) double approximationFactorForProcessMapping(const X& hypergraph, const Context& context)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(OBJECTIVE_1)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(OBJECTIVE_2)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(CONTRIBUTION)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(IS_BALANCED)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(IMBALANCE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(APPROX_FACTOR)
} // namespace mt_kahypar::metrics
