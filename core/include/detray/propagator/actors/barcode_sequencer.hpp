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

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>

namespace detray {

struct barcode_sequencer : actor {

    struct state {

        using sequence_t = vecmem::device_vector<detray::geometry::barcode>;
        sequence_t _sequence;

        /// Constructor with the vector of track states
        DETRAY_HOST_DEVICE
        explicit state(sequence_t seq) : _sequence(seq) {}
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& actor_state,
                                       propagator_state_t& propagation) const {

        const auto& navigation = propagation._navigation;

        if (!(navigation.is_on_sensitive() ||
              navigation.encountered_sf_material())) {
            return;
        }

        const auto& bcd = navigation.current().sf_desc.barcode();
        assert(!bcd.is_invalid());

        actor_state._sequence.push_back(bcd);

        return;
    }
};

}  // namespace detray
