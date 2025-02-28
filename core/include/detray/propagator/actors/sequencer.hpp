/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/barcode.hpp"
#include "detray/propagator/base_actor.hpp"

namespace detray {

struct barcode_sequencer : actor {

    struct state {

        using sequence_t = vecmem::device_vector<detray::geometry::barcode>;
        sequence_t _sequence;
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& actor_state,
                                       propagator_state_t& propagation) const {

        const auto& navigation = propagation._navigation;

        if (navigation.is_on_surface()) {
            actor_state._sequence.push_back(navigation.current());
        }

        return;
    }
};

}  // namespace detray
