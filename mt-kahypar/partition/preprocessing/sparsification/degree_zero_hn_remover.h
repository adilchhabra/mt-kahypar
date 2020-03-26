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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
template<typename TypeTraits>
class DegreeZeroHypernodeRemoverT {

  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;

 public:
  DegreeZeroHypernodeRemoverT(const Context& context) :
    _context(context),
    _removed_hns(),
    _mapping() { }

  DegreeZeroHypernodeRemoverT(const DegreeZeroHypernodeRemoverT&) = delete;
  DegreeZeroHypernodeRemoverT & operator= (const DegreeZeroHypernodeRemoverT &) = delete;

  DegreeZeroHypernodeRemoverT(DegreeZeroHypernodeRemoverT&&) = delete;
  DegreeZeroHypernodeRemoverT & operator= (DegreeZeroHypernodeRemoverT &&) = delete;

  // ! Contracts degree-zero vertices to degree-zero supervertices
  // ! We contract sets of degree-zero vertices such that the weight of
  // ! each supervertex is less than or equal than the maximum allowed
  // ! node weight for a vertex during coarsening.
  HypernodeID contractDegreeZeroHypernodes(HyperGraph& hypergraph) {
    _mapping.assign(hypergraph.initialNumNodes(), kInvalidHypernode);
    HypernodeID current_num_nodes = hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes();
    HypernodeID num_removed_degree_zero_hypernodes = 0;
    HypernodeID last_degree_zero_representative = kInvalidHypernode;
    HypernodeWeight last_degree_zero_weight = 0;
    for (const HypernodeID& hn : hypergraph.nodes()) {
      if ( current_num_nodes <= _context.coarsening.contraction_limit ) {
        break;
      }

      if ( hypergraph.nodeDegree(hn) == 0 ) {
        bool was_removed = false;
        if ( last_degree_zero_representative != kInvalidHypernode ) {
          const HypernodeWeight weight = hypergraph.nodeWeight(hn);
          if ( last_degree_zero_weight + weight <= _context.coarsening.max_allowed_node_weight ) {
            // Remove vertex and aggregate its weight in its represenative supervertex
            ++num_removed_degree_zero_hypernodes;
            --current_num_nodes;
            hypergraph.removeHypernode(hn);
            _removed_hns.push_back(hn);
            was_removed = true;
            _mapping[hypergraph.originalNodeID(hn)] = last_degree_zero_representative;
            last_degree_zero_weight += weight;
            hypergraph.setNodeWeight(last_degree_zero_representative, last_degree_zero_weight);
          }
        }

        if ( !was_removed ) {
          last_degree_zero_representative = hn;
          last_degree_zero_weight = hypergraph.nodeWeight(hn);
        }
      }
    }
    return num_removed_degree_zero_hypernodes;
  }

  // ! Restore degree-zero vertices
  // ! Each removed degree-zero vertex is assigned to the block of its supervertex.
  void restoreDegreeZeroHypernodes(PartitionedHyperGraph& hypergraph) {
    for ( const HypernodeID& hn : _removed_hns ) {
      const HypernodeID original_id = hypergraph.originalNodeID(hn);
      ASSERT(original_id < _mapping.size());
      const HypernodeID representative = _mapping[original_id];
      ASSERT(representative != kInvalidHypernode);
      // Restore degree-zero vertex and assign it to the block
      // of its supervertex
      hypergraph.enableHypernode(hn);
      hypergraph.setNodeWeight(representative,
        hypergraph.nodeWeight(representative) - hypergraph.nodeWeight(hn));
      hypergraph.setNodePart(hn, hypergraph.partID(representative));
    }
  }

  void assignAllDegreeZeroHypernodesToSameCommunity(HyperGraph& hypergraph, ds::Clustering& clustering) {
    ASSERT(hypergraph.initialNumNodes() <= clustering.size());
    PartitionID community_id = kInvalidPartition;
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      if ( hypergraph.nodeDegree(hn) == 0 ) {
        if ( community_id != kInvalidPartition ) {
          clustering[hypergraph.originalNodeID(hn)] = community_id;
        } else {
          community_id = clustering[hypergraph.originalNodeID(hn)];
        }
      }
    }
  }

 private:
  const Context& _context;
  parallel::scalable_vector<HypernodeID> _removed_hns;
  parallel::scalable_vector<HypernodeID> _mapping;
};

using DegreeZeroHypernodeRemover = DegreeZeroHypernodeRemoverT<GlobalTypeTraits>;

}  // namespace mt_kahypar