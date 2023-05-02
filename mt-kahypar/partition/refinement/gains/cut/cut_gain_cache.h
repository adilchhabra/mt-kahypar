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

#include "kahypar/meta/policy_registry.h"

#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/partition/refinement/gains/cut/cut_rollback.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/macros.h"

namespace mt_kahypar {

// Forward
class DeltaCutGainCache;

/**
 * The gain cache stores the gain values for all possible node moves for the cut metric.
 *
 * For a weighted hypergraph H = (V,E,c,w), the cut metric is defined as follows
 * connectivity(H) := \sum_{e \in cut(E)} w(e).
 *
 * The gain of moving a node u from its current block V_i to a target block V_j can be expressed as follows
 * g(u, V_j) := g(u, V_j) := w({ e \in I(u) | pin_count(e, V_j) = |e| - 1 }) - w({ e \in I(u) | pin_count(e, V_i) = |e| }).
 * Moving node u from V_i to V_j, removes all nets e \in I(u) from the cut where pin_cout(e, V_j) = |e| - 1,
 * but makes it a cut hyperedge if pin_count(e, V_i) = |e|.
 *
 * We call the first term in the equation the benefit term b(u, V_j) and the second the penalty term p(u).
 * Our gain cache stores and maintains these entries for each node and block.
 * Thus, the gain cache stores k + 1 entries per node.
*/
class CutGainCache final : public kahypar::meta::PolicyBase {

  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = ID(100000);

 public:
  static constexpr GainPolicy TYPE = GainPolicy::cut;
  using DeltaGainCache = DeltaCutGainCache;
  using Rollback = CutRollback;

  CutGainCache() :
    _is_initialized(false),
    _k(kInvalidPartition),
    _gain_cache() { }

  CutGainCache(const CutGainCache&) = delete;
  CutGainCache & operator= (const CutGainCache &) = delete;

  CutGainCache(CutGainCache&& other) = default;
  CutGainCache & operator= (CutGainCache&& other) = default;

  // ####################### Initialization #######################

  bool isInitialized() const {
    return _is_initialized;
  }

  void reset(const bool run_parallel = true) {
    unused(run_parallel);
    _is_initialized = false;
  }

  size_t size() const {
    return _gain_cache.size();
  }

  // ! Initializes all gain cache entries
  template<typename PartitionedHypergraph>
  void initializeGainCache(const PartitionedHypergraph& partitioned_hg);

  // ####################### Gain Computation #######################

  // ! Returns the penalty term of node u.
  // ! More formally, p(u) := (w(I(u)) - w({ e \in I(u) | pin_count(e, V_i) = |e| })
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight penaltyTerm(const HypernodeID u,
                              const PartitionID /* only relevant for graphs */) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return _gain_cache[penalty_index(u)].load(std::memory_order_relaxed);
  }

  // ! Recomputes the penalty term entry in the gain cache
  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void recomputePenaltyTermEntry(const PartitionedHypergraph& partitioned_hg,
                                 const HypernodeID u) {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    _gain_cache[penalty_index(u)].store(recomputePenaltyTerm(
      partitioned_hg, u), std::memory_order_relaxed);
  }

  // ! Returns the benefit term for moving node u to block to.
  // ! More formally, b(u, V_j) := w({ e \in I(u) | pin_count(e, V_j) = |e| - 1 })
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight benefitTerm(const HypernodeID u, const PartitionID to) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return _gain_cache[benefit_index(u, to)].load(std::memory_order_relaxed);
  }

  // ! Returns the gain of moving node u from its current block to a target block V_j.
  // ! More formally, g(u, V_j) := b(u, V_j) - p(u).
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight gain(const HypernodeID u,
                       const PartitionID, /* only relevant for graphs */
                       const PartitionID to ) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return benefitTerm(u, to) - penaltyTerm(u, kInvalidPartition);
  }

  // ####################### Delta Gain Update #######################

  // ! This functions implements the delta gain updates for the cut metric.
  // ! When moving a node from its current block from to a target block to, we iterate
  // ! over its incident hyperedges and update their pin count values. After each pin count
  // ! update, we call this function to update the gain cache to changes associated with
  // ! corresponding hyperedge.
  template<typename PartitionedHypergraph>
  void deltaGainUpdate(const PartitionedHypergraph& partitioned_hg,
                       const HyperedgeID he,
                       const HyperedgeWeight we,
                       const PartitionID from,
                       const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to,
                       const HypernodeID pin_count_in_to_part_after);

  // ! Computes the improvement for cut metric associated with a hyperedge.
  static HyperedgeWeight delta(const HyperedgeID,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
    return edge_size > 1 ? (pin_count_in_from_part_after == edge_size - 1) * edge_weight -
      ( pin_count_in_to_part_after == edge_size ) * edge_weight : 0;
  }

  // ####################### Uncontraction #######################

  // ! This function implements the gain cache update after an uncontraction that restores node v in
  // ! hyperedge he. After the uncontraction node u and v are contained in hyperedge he.
  template<typename PartitionedHypergraph>
  void uncontractUpdateAfterRestore(const PartitionedHypergraph& partitioned_hg,
                                    const HypernodeID u,
                                    const HypernodeID v,
                                    const HyperedgeID he,
                                    const HypernodeID pin_count_in_part_after);

  // ! This function implements the gain cache update after an uncontraction that replaces u with v in
  // ! hyperedge he. After the uncontraction only node v is contained in hyperedge he.
  template<typename PartitionedHypergraph>
  void uncontractUpdateAfterReplacement(const PartitionedHypergraph& partitioned_hg,
                                        const HypernodeID u,
                                        const HypernodeID v,
                                        const HyperedgeID he);

  // ! This function is called after restoring a single-pin hyperedge. The function assumes that
  // ! u is the only pin of the corresponding hyperedge, while block_of_u is its corresponding block ID.
  void restoreSinglePinHyperedge(const HypernodeID,
                                 const PartitionID,
                                 const HyperedgeWeight) {
    // Do nothing here
  }

  // ####################### Only for Testing #######################

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight recomputePenaltyTerm(const PartitionedHypergraph& partitioned_hg,
                                       const HypernodeID u) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    const PartitionID block_of_u = partitioned_hg.partID(u);
    HyperedgeWeight penalty = 0;
    for (HyperedgeID e : partitioned_hg.incidentEdges(u)) {
      const HypernodeID edge_size = partitioned_hg.edgeSize(e);
      if ( partitioned_hg.pinCountInPart(e, block_of_u) == edge_size && edge_size > 1 ) {
        penalty += partitioned_hg.edgeWeight(e);
      }
    }
    return penalty;
  }

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight recomputeBenefitTerm(const PartitionedHypergraph& partitioned_hg,
                                       const HypernodeID u,
                                       const PartitionID to) const {
    HyperedgeWeight benefit = 0;
    for (HyperedgeID e : partitioned_hg.incidentEdges(u)) {
      const HypernodeID edge_size = partitioned_hg.edgeSize(e);
      if (partitioned_hg.pinCountInPart(e, to) == edge_size - 1 && edge_size > 1) {
        benefit += partitioned_hg.edgeWeight(e);
      }
    }
    return benefit;
  }

 private:
  friend class DeltaCutGainCache;

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t penalty_index(const HypernodeID u) const {
    return size_t(u) * ( _k + 1 );
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t benefit_index(const HypernodeID u, const PartitionID p) const {
    return size_t(u) * ( _k + 1 )  + p + 1;
  }

  // ! Allocates the memory required to store the gain cache
  void allocateGainTable(const HypernodeID num_nodes,
                         const PartitionID k) {
    if (_gain_cache.size() == 0) {
      ASSERT(_k == kInvalidPartition);
      _k = k;
      _gain_cache.resize(
        "Refinement", "gain_cache", num_nodes * size_t(_k + 1), true);
    }
  }

  // ! Initializes the benefit and penalty terms for a node u
  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void initializeGainCacheEntryForNode(const PartitionedHypergraph& partitioned_hg,
                                       const HypernodeID u,
                                       vec<Gain>& benefit_aggregator);

  bool nodeGainAssertions(const HypernodeID u, const PartitionID p) const {
    if ( p == kInvalidPartition || p >= _k ) {
      LOG << "Invalid block ID (Node" << u << "is part of block" << p
          << ", but valid block IDs must be in the range [ 0," << _k << "])";
      return false;
    }
    if ( benefit_index(u, p) >= _gain_cache.size() ) {
      LOG << "Access to gain cache would result in an out-of-bounds access ("
          << "Benefit Index =" << benefit_index(u, p)
          << ", Gain Cache Size =" << _gain_cache.size() << ")";
      return false;
    }
    return true;
  }


  // ! Indicate whether or not the gain cache is initialized
  bool _is_initialized;

  // ! Number of blocks
  PartitionID _k;

  // ! Array of size |V| * (k + 1), which stores the benefit and penalty terms of each node.
  ds::Array< CAtomic<HyperedgeWeight> > _gain_cache;
};

/**
 * In our FM algorithm, the different local searches perform nodes moves locally not visible for other
 * threads. The delta gain cache stores these local changes relative to the shared
 * gain cache. For example, the penalty term can be computed as follows
 * p'(u) := p(u) + Δp(u)
 * where p(u) is the penalty term stored in the shared gain cache and Δp(u) is the penalty term stored in
 * the delta gain cache after performing some moves locally. To maintain Δp(u) and Δb(u,V_j), we use a hash
 * table that only stores entries affected by a gain cache update.
*/
class DeltaCutGainCache {

 public:
  DeltaCutGainCache(const CutGainCache& gain_cache) :
    _gain_cache(gain_cache),
    _gain_cache_delta() { }

  // ####################### Initialize & Reset #######################

  void initialize(const size_t size) {
    _gain_cache_delta.initialize(size);
  }

  void clear() {
    _gain_cache_delta.clear();
  }

  void dropMemory() {
    _gain_cache_delta.freeInternalData();
  }

  size_t size_in_bytes() const {
    return _gain_cache_delta.size_in_bytes();
  }

  // ####################### Gain Computation #######################

  // ! Returns the penalty term of node u.
  // ! More formally, p(u) := (w(I(u)) - w({ e \in I(u) | pin_count(e, V_i) = |e| })
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight penaltyTerm(const HypernodeID u,
                              const PartitionID from) const {
    const HyperedgeWeight* penalty_delta =
      _gain_cache_delta.get_if_contained(_gain_cache.penalty_index(u));
    return _gain_cache.penaltyTerm(u, from) + ( penalty_delta ? *penalty_delta : 0 );
  }

  // ! Returns the benefit term for moving node u to block to.
  // ! More formally, b(u, V_j) := w({ e \in I(u) | pin_count(e, V_j) = |e| - 1 })
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight benefitTerm(const HypernodeID u, const PartitionID to) const {
    ASSERT(to != kInvalidPartition && to < _gain_cache._k);
    const HyperedgeWeight* benefit_delta =
      _gain_cache_delta.get_if_contained(_gain_cache.benefit_index(u, to));
    return _gain_cache.benefitTerm(u, to) + ( benefit_delta ? *benefit_delta : 0 );
  }

  // ! Returns the gain of moving node u from its current block to a target block V_j.
  // ! More formally, g(u, V_j) := b(u, V_j) - p(u).
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight gain(const HypernodeID u,
                       const PartitionID from,
                       const PartitionID to ) const {
    return benefitTerm(u, to) - penaltyTerm(u, from);
  }

 // ####################### Delta Gain Update #######################

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void deltaGainUpdate(const PartitionedHypergraph& partitioned_hg,
                       const HyperedgeID he,
                       const HyperedgeWeight we,
                       const PartitionID from,
                       const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to,
                       const HypernodeID pin_count_in_to_part_after) {
    const HypernodeID edge_size = partitioned_hg.edgeSize(he);
    if ( edge_size > 1 ) {
      if ( pin_count_in_from_part_after == edge_size - 1 ) {
        for ( const HypernodeID& u : partitioned_hg.pins(he) ) {
          _gain_cache_delta[_gain_cache.penalty_index(u)] -= we;
          _gain_cache_delta[_gain_cache.benefit_index(u, from)] += we;
        }
      } else if ( pin_count_in_from_part_after == edge_size - 2 ) {
        for ( const HypernodeID& u : partitioned_hg.pins(he) ) {
          _gain_cache_delta[_gain_cache.benefit_index(u, from)] -= we;
        }
      }

      if ( pin_count_in_to_part_after == edge_size ) {
        for ( const HypernodeID& u : partitioned_hg.pins(he) ) {
          _gain_cache_delta[_gain_cache.penalty_index(u)] += we;
          _gain_cache_delta[_gain_cache.benefit_index(u, to)] -= we;
        }
      } else if ( pin_count_in_to_part_after == edge_size - 1 ) {
        for ( const HypernodeID& u : partitioned_hg.pins(he) ) {
          _gain_cache_delta[_gain_cache.benefit_index(u, to)] += we;
        }
      }
    }
  }

 // ####################### Miscellaneous #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
    utils::MemoryTreeNode* gain_cache_delta_node = parent->addChild("Delta Gain Cache");
    gain_cache_delta_node->updateSize(size_in_bytes());
  }

 private:
  const CutGainCache& _gain_cache;

  // ! Stores the delta of each locally touched gain cache entry
  // ! relative to the gain cache in '_phg'
  ds::DynamicFlatMap<size_t, HyperedgeWeight> _gain_cache_delta;
};

}  // namespace mt_kahypar