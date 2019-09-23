/*******************************************************************************
 * This file is part of KaHyPar.
 *
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

#include <iostream>
#include <string>
  
#include "kahypar/macros.h"

namespace mt_kahypar {

enum class Type : int8_t {
  Unweighted = 0,
  EdgeWeights = 1,
  NodeWeights = 10,
  EdgeAndNodeWeights = 11,
};

enum class CommunityAssignmentObjective {
  vertex_objective,
  pin_objective,
  undefined
};

enum class CommunityAssignmentStrategy {
  bin_packing,
  undefined
};


std::ostream& operator<< (std::ostream& os, const Type& type) {
  switch (type) {
    case Type::Unweighted: return os << "unweighted";
    case Type::EdgeWeights: return os << "edge_weights";
    case Type::NodeWeights: return os << "node_weights";
    case Type::EdgeAndNodeWeights: return os << "edge_and_node_weights";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(type);
}

std::ostream& operator<< (std::ostream& os, const CommunityAssignmentObjective& objective) {
  switch (objective) {
    case CommunityAssignmentObjective::vertex_objective: return os << "vertex_objective";
    case CommunityAssignmentObjective::pin_objective: return os << "pin_objective";
    case CommunityAssignmentObjective::undefined: return os << "undefined";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(objective);
}

std::ostream& operator<< (std::ostream& os, const CommunityAssignmentStrategy& strategy) {
  switch (strategy) {
    case CommunityAssignmentStrategy::bin_packing: return os << "bin_packing";
    case CommunityAssignmentStrategy::undefined: return os << "undefined";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(strategy);
}

static CommunityAssignmentObjective communityAssignmentObjectiveFromString(const std::string& objective) {
  if (objective == "vertex_objective") {
    return CommunityAssignmentObjective::vertex_objective;
  } else if (objective == "pin_objective") {
    return CommunityAssignmentObjective::pin_objective;
  }
  LOG << "No valid community assignment objective.";
  exit(0);
  return CommunityAssignmentObjective::undefined;
}

static CommunityAssignmentStrategy communityAssignmentStrategyFromString(const std::string& objective) {
  if (objective == "bin_packing") {
    return CommunityAssignmentStrategy::bin_packing;
  }
  LOG << "No valid community assignment strategy.";
  exit(0);
  return CommunityAssignmentStrategy::undefined;
}

} // namesapce mt_kahypar