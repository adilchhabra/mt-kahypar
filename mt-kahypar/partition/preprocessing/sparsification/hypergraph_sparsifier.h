/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <vector>

#include "kahypar/utils/math.h"
#include "kahypar/utils/hash_vector.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/datastructures/streaming_map.h"
#include "mt-kahypar/datastructures/sparsifier_hypergraph.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/preprocessing/sparsification/i_hypergraph_sparsifier.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
template<typename TypeTraits,
         typename SimiliarNetCombiner>
class HypergraphSparsifierT : public IHypergraphSparsifierT<TypeTraits> {

  using Base = IHypergraphSparsifierT<TypeTraits>;
  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using HyperGraphFactory = typename TypeTraits::HyperGraphFactory;
  using TBB = typename TypeTraits::TBB;
  using SparsifierHypergraph = ds::SparsifierHypergraph<HyperGraph, HyperGraphFactory, TBB>;

  using HashFunc = kahypar::math::MurmurHash<HypernodeID>;
  using HashValue = typename HashFunc::HashValue;
  using HashFuncVector = kahypar::HashFuncVector<HashFunc>;

  struct Footprint {
    explicit Footprint() :
      footprint(),
      he(kInvalidHyperedge) { }

    parallel::scalable_vector<HashValue> footprint;
    HyperedgeID he;

    bool operator==(const Footprint& other) {
      ASSERT(footprint.size() == other.footprint.size());
      for ( size_t i = 0; i < footprint.size(); ++i ) {
        if ( footprint[i] != other.footprint[i] ) {
          return false;
        }
      }
      return true;
    }

    bool operator<(const Footprint& other) {
      ASSERT(footprint.size() == other.footprint.size());
      for ( size_t i = 0; i < footprint.size(); ++i ) {
        if ( footprint[i] < other.footprint[i] ) {
          return true;
        } else if ( footprint[i] > other.footprint[i] ) {
          return false;
        }
      }
      return he < other.he;
    }

  };

  using FootprintMap = parallel::scalable_vector<parallel::scalable_vector<Footprint>>;

  static constexpr bool enable_heavy_assert = false;

 public:
  HypergraphSparsifierT(const Context& context,
                        const TaskGroupID task_group_id) :
    Base(),
    _context(context),
    _task_group_id(task_group_id),
    _sparsified_hg(),
    _sparsified_partitioned_hg(),
    _mapping() { }

  HypergraphSparsifierT(const HypergraphSparsifierT&) = delete;
  HypergraphSparsifierT & operator= (const HypergraphSparsifierT &) = delete;

  HypergraphSparsifierT(HypergraphSparsifierT&&) = delete;
  HypergraphSparsifierT & operator= (HypergraphSparsifierT &&) = delete;


 private:

  // ####################### Sparsification Functions #######################

  HyperGraph& sparsifiedHypergraphImpl() override final {
    return _sparsified_hg;
  }

  PartitionedHyperGraph& sparsifiedPartitionedHypergraphImpl() override final {
    return _sparsified_partitioned_hg;
  }

  void sparsifyImpl(const HyperGraph& hypergraph) override final {
    ASSERT(_context.useSparsification());
    SparsifierHypergraph sparsified_hypergraph(hypergraph, _task_group_id);

    // #################### STAGE 1 ####################
    // Heavy Hyperedge Removal
    // If the weight of all pins of a hyperedge is greater than a
    // certain threshold, we remove them from the hypergraph
    if ( _context.sparsification.use_heavy_net_removal ) {
      utils::Timer::instance().start_timer("heavy_hyperedge_removal", "Heavy HE Removal");
      heavyHyperedgeRemovalSparsification(sparsified_hypergraph);
      utils::Timer::instance().stop_timer("heavy_hyperedge_removal");
    }

    // #################### STAGE 2 ####################
    if ( _context.sparsification.use_similiar_net_removal ) {
      utils::Timer::instance().start_timer("similiar_hyperedge_removal", "Similiar HE Removal");
      similiarHyperedgeRemoval(hypergraph, sparsified_hypergraph);
      utils::Timer::instance().stop_timer("similiar_hyperedge_removal");
    }

    // #################### STAGE 3 ####################
    // Perform Degree-Zero Contractions
    // Degree-Zero hypernodes are contracted to supervertices such that
    // each supervertex has a weight smaller than the maximum allowed
    // node weight.
    if ( _context.sparsification.use_degree_zero_contractions ) {
      utils::Timer::instance().start_timer("degree_zero_contraction", "Degree-Zero Contractions");
      degreeZeroSparsification(sparsified_hypergraph);
      utils::Timer::instance().stop_timer("degree_zero_contraction");
    }

    // #################### STAGE 4 ####################
    // Construct sparsified hypergraph
    utils::Timer::instance().start_timer("construct_sparsified_hypergraph", "Construct Sparsified HG");
    _sparsified_hg = sparsified_hypergraph.sparsify();
    _mapping = sparsified_hypergraph.getMapping();
    _sparsified_partitioned_hg = PartitionedHyperGraph(
      _context.partition.k, _task_group_id, _sparsified_hg);
    utils::Timer::instance().stop_timer("construct_sparsified_hypergraph");
  }

  void undoSparsificationImpl(PartitionedHyperGraph& hypergraph) override final {
    hypergraph.doParallelForAllNodes(_task_group_id, [&](const HypernodeID hn) {
      const HypernodeID original_id = hypergraph.originalNodeID(hn);
      ASSERT(original_id < _mapping.size());
      const HypernodeID original_sparsified_id = _mapping[original_id];
      const HypernodeID sparsified_hn = _sparsified_partitioned_hg.globalNodeID(original_sparsified_id);
      ASSERT(_sparsified_partitioned_hg.nodeIsEnabled(sparsified_hn));
      hypergraph.setNodePart(hn, _sparsified_partitioned_hg.partID(sparsified_hn));
    });
    hypergraph.initializeNumCutHyperedges(_task_group_id);
  }

 private:
  // ! Similiar to contractDegreeZeroHypernodes, but contraction is not applied
  // ! directly to hypergraph. Instead a mapping is computed that maps each vertex
  // ! of the original hypergraph to its supervertex and the weight of each
  // ! supervertex is aggregated in the hypernode weight vector.
  void degreeZeroSparsification(SparsifierHypergraph& hypergraph) {
    HypernodeID current_num_nodes = hypergraph.numNodes() - hypergraph.numRemovedNodes();
    HypernodeID degree_zero_supervertex = kInvalidHypernode;
    for (HypernodeID hn = 0; hn < hypergraph.numNodes(); ++hn) {
      if ( current_num_nodes <= _context.coarsening.contraction_limit ) {
        break;
      }

      if ( hypergraph.nodeDegree(hn) == 0 ) {
        bool was_removed = false;
        if ( degree_zero_supervertex != kInvalidHypernode ) {
          if ( hypergraph.nodeWeight(degree_zero_supervertex) +
               hypergraph.nodeWeight(hn) <=
               _context.coarsening.max_allowed_node_weight ) {
            // Remove vertex and aggregate its weight in its represenative supervertex
            hypergraph.contract(degree_zero_supervertex, hn);
            --current_num_nodes;
            was_removed = true;
          }
        }

        if ( !was_removed ) {
          degree_zero_supervertex = hn;
        }
      }
    }
  }

  // ! Removes hyperedges where the weight of all pins is greater
  // ! than a certain threshold. The threshold is specified in
  // ! '_context.initial_partitioning.max_hyperedge_pin_weight'
  void heavyHyperedgeRemovalSparsification(SparsifierHypergraph& hypergraph) {
    tbb::parallel_for(ID(0), hypergraph.numEdges(), [&](const HyperedgeID& e) {
      if ( hypergraph.edgeIsEnabled(e) ) {
        HypernodeWeight pin_weight = 0;
        for ( const HypernodeID& pin : hypergraph.pins(e) ) {
          pin_weight += hypergraph.nodeWeight(pin);
        }
        // Hyperedge will be include in sparsified hypergraph if its weight of
        // all pins is less than a predefined upper bound
        if ( pin_weight >= _context.sparsification.max_hyperedge_pin_weight ) {
          hypergraph.remove(e);
        }
      }
    });
  }

  void similiarHyperedgeRemoval(const HyperGraph& original_hg, SparsifierHypergraph& hypergraph) {
    HashFuncVector hash_functions(_context.sparsification.min_hash_footprint_size,
      utils::Randomize::instance().getRandomInt(0, 1000, sched_getcpu()));

    ds::StreamingMap<HashValue, Footprint> hash_buckets;
    tbb::parallel_for(ID(0), hypergraph.numEdges(), [&](const HyperedgeID he) {
      if ( hypergraph.edgeIsEnabled(he) ) {
        Footprint he_footprint;
        he_footprint.footprint = {};
        he_footprint.he = he;
        for ( size_t i = 0; i < hash_functions.getHashNum(); ++i ) {
          he_footprint.footprint.push_back(minHash(hash_functions[i], hypergraph.pins(he)));
        }
        hash_buckets.stream(combineHash(he_footprint), std::move(he_footprint));
      }
    });

    FootprintMap footprint_map(hash_buckets.size());
    hash_buckets.copy(footprint_map, [&](const HashValue key) {
      return key % hash_buckets.size();
    });

    tbb::parallel_for(0UL, footprint_map.size(), [&](const size_t bucket) {
      parallel::scalable_vector<Footprint>& footprint_bucket = footprint_map[bucket];
      if ( footprint_bucket.size() > 0 ) {
        std::sort(footprint_bucket.begin(), footprint_bucket.end());

        for ( size_t i = 0; i < footprint_bucket.size(); ++i ) {
          Footprint& representative = footprint_bucket[i];
          if ( representative.he != kInvalidHyperedge ) {
            parallel::scalable_vector<HypernodeID> rep_he = hypergraph.pins(representative.he);
            HyperedgeWeight rep_weight = hypergraph.edgeWeight(representative.he);
            bool exist_similiar_hes = false;
            for ( size_t j = i + 1; j < footprint_bucket.size(); ++j ) {
              Footprint& similiar_footprint = footprint_bucket[j];
              if ( similiar_footprint.he != kInvalidHyperedge ) {
                if ( representative == similiar_footprint ) {
                  const double jaccard_index = jaccard(
                    hypergraph.pins(representative.he), hypergraph.pins(similiar_footprint.he));
                  if ( jaccard_index >= _context.sparsification.jaccard_threshold ) {
                    rep_he = SimiliarNetCombiner::combine(original_hg, rep_he, hypergraph.pins(similiar_footprint.he));
                    rep_weight += hypergraph.edgeWeight(similiar_footprint.he);
                    hypergraph.remove(similiar_footprint.he);
                    similiar_footprint.he = kInvalidHyperedge;
                    exist_similiar_hes = true;
                  }
                } else {
                  break;
                }
              }
            }

            if ( exist_similiar_hes ) {
              hypergraph.replace(representative.he, std::move(rep_he));
              hypergraph.setEdgeWeight(representative.he, rep_weight);
            }
          }
        }
      }
    });
  }

  HashValue minHash(const HashFunc& hash_function,
                    const parallel::scalable_vector<HypernodeID>& hyperedge ) {
    HashValue hash_value = std::numeric_limits<HashValue>::max();
    for ( const HypernodeID& pin : hyperedge ) {
      hash_value = std::min(hash_value, hash_function(pin));
    }
    return hash_value;
  }

  HashValue combineHash(const Footprint& footprint) {
    HashValue hash_value = kEdgeHashSeed;
    for ( const HashValue& value : footprint.footprint ) {
      hash_value ^= value;
    }
    return hash_value;
  }

  double jaccard(const parallel::scalable_vector<HypernodeID>& lhs,
                 const parallel::scalable_vector<HypernodeID>& rhs) {
    const size_t min_size = std::min(lhs.size(), rhs.size());
    const size_t max_size = std::max(lhs.size(), rhs.size());
    if ( static_cast<double>(min_size) / static_cast<double>(max_size) <
         _context.sparsification.jaccard_threshold ) {
      return 0.0;
    }

    size_t intersection_size = 0;
    size_t i = 0;
    size_t j = 0;
    while ( i < lhs.size() && j < rhs.size() ) {
      if ( lhs[i] == rhs[j] ) {
        ++intersection_size;
        ++i;
        ++j;
      } else if ( lhs[i] < rhs[j] ) {
        ++i;
      } else {
        ++j;
      }
    }
    const size_t union_size = lhs.size() + rhs.size() - intersection_size;
    return static_cast<double>(intersection_size) /
      static_cast<double>(union_size);
  }

  const Context& _context;
  const TaskGroupID _task_group_id;

  HyperGraph _sparsified_hg;
  PartitionedHyperGraph _sparsified_partitioned_hg;
  parallel::scalable_vector<HypernodeID> _mapping;
};

template<typename SimiliarNetCombiner>
using HypergraphSparsifier = HypergraphSparsifierT<GlobalTypeTraits, SimiliarNetCombiner>;

}  // namespace mt_kahypar