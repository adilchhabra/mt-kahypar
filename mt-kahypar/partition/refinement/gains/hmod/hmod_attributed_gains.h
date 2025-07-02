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

#include "mt-kahypar/datastructures/hypergraph_common.h"

namespace mt_kahypar {

/**
 * After moving a node, we perform a synchronized update of the pin count values
 * for each incident hyperedge of the node based on which we then compute an
 * attributed gain value.
 */
struct HModAttributedGains {
  static Gain gain(const SynchronizedEdgeUpdate &sync_update) {
    const HypernodeID d = static_cast<HypernodeID>(sync_update.edge_size);
    if (d <= 1)
      return 0;

    /* ---------- 1. ΔEC for this edge ------------------------------------ */
    double deltaEC = 0.0;
    const HypernodeID old_from =
        static_cast<HypernodeID>(sync_update.pin_count_in_from_part_after + 1); // before move
    if (old_from == d) {
      deltaEC -= sync_update.edge_weight;
    } // edge *leaves* from
    if (sync_update.pin_count_in_to_part_after == d) {
      deltaEC += sync_update.edge_weight;
    }                                              // edge *enters* to
    //deltaEC /= static_cast<double>(sync_update.m); // normalise

    // --- degree-tax part -----------------------------------------------
    // DT change is a *global* quantity.  We attribute it evenly
    // across all incident edges of hn.
    //
    //   ΔDT = DT(vol_to+str) - DT(vol_to)   - ( DT(vol_from) - DT(vol_from-str)
    //   )
    //
    // Those four volumes are already available in the synchronized update.
    //
    //   Here we distribute that total change uniformly:
    auto DT = [&](double V) {
      return V / static_cast<double>(sync_update.vol_H_top); // strict-τ closed form
    };

    const double total_deltaDT =
        DT(sync_update.vol_To) -
        DT(sync_update.vol_To - sync_update.hn_volume) -
        (DT(sync_update.vol_From + sync_update.hn_volume) -
         DT(sync_update.vol_From));

    const double deltaDT_edge =
        total_deltaDT / sync_update.hn_degree; // equal share

    /* ---------- 3. final ΔQ for this edge ------------------------------- */
    const double resolution = 1.0; // γ = resolution here
    const double deltaQ = deltaEC - resolution * deltaDT_edge;

    //return static_cast<HyperedgeWeight>(std::llround(-deltaQ * 100));
    return -1.0 * deltaQ;

  }
};
} // namespace mt_kahypar
