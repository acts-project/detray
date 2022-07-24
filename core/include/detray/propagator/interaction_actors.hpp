/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/materials/interaction_kernel.hpp"
#include "detray/propagator/base_actor.hpp"

namespace detray {

struct pointwise_material_interactor : actor {

    struct state {};

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state &interactor_state,
                                       propagator_state_t &prop_state) const {}
};

}  // namespace detray