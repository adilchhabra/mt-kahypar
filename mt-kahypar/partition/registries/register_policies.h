/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2018 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "kahypar/meta/policy_registry.h"
#include "kahypar/meta/registrar.h"

#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/partition/preprocessing/policies/community_assignment_objective.h"

#define REGISTER_POLICY(policy, id, policy_class)                                  \
  static meta::Registrar<meta::PolicyRegistry<policy> > register_ ## policy_class( \
    id, new policy_class())

namespace mt_kahypar {

// //////////////////////////////////////////////////////////////////////////////
//                        Community Assignment Strategy
// //////////////////////////////////////////////////////////////////////////////
REGISTER_POLICY(CommunityAssignmentObjective, CommunityAssignmentObjective::vertex_objective,
                VertexObjectivePolicy);
REGISTER_POLICY(CommunityAssignmentObjective, CommunityAssignmentObjective::pin_objective,
                PinObjectivePolicy);

} // namespace mt_kahypar