/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2024
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

#include "gmock/gmock.h"

#include <cmath>
#include <limits>
#include <vector>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/gains/aon_hypermodularity/aon_hypermodularity_gain_computation.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

namespace {
using TypeTraits = StaticHypergraphTypeTraits;
using Hypergraph = typename TypeTraits::Hypergraph;
using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
using HypergraphFactory = typename Hypergraph::Factory;

auto make_context(const PartitionID k) {
  Context ctx;
  ctx.partition.k = k;
  ctx.partition.objective = Objective::aon_hypermodularity;
  ctx.partition.max_part_weights.assign(k, std::numeric_limits<HypernodeWeight>::max());
  return ctx;
}

void assign_partition_ids(Hypergraph& hg,
                          PartitionedHypergraph& phg,
                          const std::vector<PartitionID>& part_ids) {
  for (HypernodeID u = 0; u < part_ids.size(); ++u) {
    hg.setCommunityID(u, part_ids[u]);
  }
  hg.computeAONParameters();
  for (HypernodeID u = 0; u < part_ids.size(); ++u) {
    phg.setNodePart(u, part_ids[u]);
  }
}

}  // namespace

TEST(AONHyperModularityTest, SyncUpdateCarriesTopLevelEdgeInfo) {
  Hypergraph hg = HypergraphFactory::construct(
    4, 5, { {0, 1, 2}, {2, 3}, {0, 3}, {0, 1}, {1, 2, 3} });
  auto context = make_context(2);
  PartitionedHypergraph phg(context.partition.k, hg, parallel_tag_t());
  const std::vector<PartitionID> part_ids = { 0, 0, 0, 1 };
  assign_partition_ids(hg, phg, part_ids);

  AONHyperModularityGainComputation gain(context, true);

  const HyperedgeID monitored_edge = 0; // {0,1,2}
  bool saw_expected = false;
  SynchronizedEdgeUpdate recorded;
  ASSERT_TRUE(phg.changeNodePart(0, 0, 1,
    [&](const SynchronizedEdgeUpdate& sync_update) {
      gain.computeDeltaForHyperedge(sync_update);
      if (sync_update.he == monitored_edge && !saw_expected) {
        recorded = sync_update;
        saw_expected = true;
      }
    }));
  ASSERT_TRUE(saw_expected);

  EXPECT_EQ(recorded.edge_size, hg.edgeSize(monitored_edge));
  EXPECT_EQ(recorded.edge_strength, hg.topLevelEdgeSize(monitored_edge));
  EXPECT_EQ(recorded.edge_weight, hg.edgeWeight(monitored_edge));
  EXPECT_DOUBLE_EQ(recorded.hn_volume, phg.nodeVolume(recorded.hn));
  EXPECT_DOUBLE_EQ(recorded.hn_strength, phg.nodeStrength(recorded.hn));
  EXPECT_EQ(recorded.beta_vec, &phg.betaVector());
  EXPECT_EQ(recorded.gamma_vec, &phg.gammaVector());
  EXPECT_EQ(recorded.max_edge_size, phg.topLevelMaxEdgeSize());
  EXPECT_NE(recorded.beta_vec, nullptr);
  EXPECT_NE(recorded.gamma_vec, nullptr);
  EXPECT_DOUBLE_EQ((*recorded.beta_vec)[recorded.edge_strength],
                   hg.beta(recorded.edge_strength));
  EXPECT_DOUBLE_EQ((*recorded.gamma_vec)[recorded.edge_strength],
                   hg.gamma(recorded.edge_strength));
  EXPECT_DOUBLE_EQ(recorded.vol_From, phg.partVolume(recorded.from));
  EXPECT_DOUBLE_EQ(recorded.vol_To, phg.partVolume(recorded.to));
}

TEST(AONHyperModularityTest, TopLevelEdgeSizeSurvivesCoarsening) {
  Hypergraph hg = HypergraphFactory::construct(
    4, 2, { {0, 1, 2}, {2, 3} });
  parallel::scalable_vector<HypernodeID> communities(hg.initialNumNodes());
  communities[0] = 0;
  communities[1] = 0; // contract node 1 with node 0
  communities[2] = 1;
  communities[3] = 2;

  Hypergraph coarse = hg.contract(communities, false);
  bool saw_shrunk_edge = false;
  for (const HyperedgeID he : coarse.edges()) {
    if (!coarse.edgeIsEnabled(he)) {
      continue;
    }
    const HypernodeID current_size = coarse.edgeSize(he);
    const HypernodeID top_level_size = coarse.topLevelEdgeSize(he);
    EXPECT_GE(top_level_size, current_size);
    if (top_level_size > current_size) {
      saw_shrunk_edge = true;
    }
  }
  EXPECT_TRUE(saw_shrunk_edge);
}

TEST(AONHyperModularityTest, MatchesReferenceBetasGammasAndObjective) {
  Hypergraph hg = HypergraphFactory::construct(
    4, 5, { {0, 1, 2}, {2, 3}, {0, 3}, {0, 1}, {1, 2, 3} });
  auto context = make_context(2);
  PartitionedHypergraph phg(context.partition.k, hg, parallel_tag_t());
  const std::vector<PartitionID> part_ids = { 0, 0, 0, 1 };
  assign_partition_ids(hg, phg, part_ids);

  ASSERT_GE(hg.maxEdgeSize(), HyperedgeID(3));
  const double expected_beta_2 = -1.194022473472768;
  const double expected_beta_3 = 0.26126475913407354;
  const double expected_gamma_2 = -0.025814814814814815;
  const double expected_gamma_3 = 0.00030717225161669567;
  EXPECT_NEAR(hg.beta(2), expected_beta_2, 1e-12);
  EXPECT_NEAR(hg.beta(3), expected_beta_3, 1e-12);
  EXPECT_NEAR(hg.gamma(2), expected_gamma_2, 1e-12);
  EXPECT_NEAR(hg.gamma(3), expected_gamma_3, 1e-12);

  const double expected_quality = -2.1267801878114625;
  const double actual_quality = metrics::quality(phg, context, false);
  EXPECT_NEAR(actual_quality, expected_quality, 1e-12);
}

}  // namespace mt_kahypar
