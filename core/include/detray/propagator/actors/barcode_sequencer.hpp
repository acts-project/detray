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
#include "detray/utils/log.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>

namespace detray {

template <typename sf_descriptor_t>
struct barcode_sequencer : actor {

    struct state {
        using surface_type = sf_descriptor_t;
        using sequence_t = vecmem::device_vector<sf_descriptor_t>;

        /// Constructor with the vector of track states
        DETRAY_HOST_DEVICE
        explicit state(sequence_t seq) : m_sequence(seq) {}

        /// Flags that the internal buffer has experienced an overflow
        DETRAY_HOST_DEVICE void set_overflow() { m_overflow = true; }

        /// @returns true if the internal buffer has reached its capacity
        DETRAY_HOST_DEVICE bool overflow() const { return m_overflow; }

        /// @returns access to the surface sequence buffer
        DETRAY_HOST_DEVICE const sequence_t& sequence() const {
            return m_sequence;
        }

        /// @returns access to the surface sequence buffer
        DETRAY_HOST_DEVICE sequence_t& sequence() { return m_sequence; }

        private:
        /// Sequence of surfaces the navigation encountered
        sequence_t m_sequence;
        /// Internal buffer overflow
        bool m_overflow = false;
    };

    template <typename propagator_state_t>
        requires std::same_as<
            sf_descriptor_t,
            typename propagator_state_t::detector_type::surface_type>
    DETRAY_HOST_DEVICE void operator()(state& actor_state,
                                       propagator_state_t& propagation) const {

        auto& navigation = propagation._navigation;

        DETRAY_VERBOSE_HOST_DEVICE(
            "Checking for next sensitive in sequence...");
        if (!(navigation.is_on_sensitive() ||
              navigation.encountered_sf_material())) {
            return;
        }

        if (actor_state.sequence().size() ==
            actor_state.sequence().capacity()) {
            DETRAY_ERROR_HOST_DEVICE("Sequence overflow!");
            actor_state.set_overflow();
            navigation.exit();
            return;
        }

        const auto& sf_desc = std::as_const(navigation).current().sf_desc;
        assert(!sf_desc.barcode().is_invalid());

        actor_state.sequence().push_back(sf_desc);
        DETRAY_VERBOSE_HOST("Added: " << sf_desc);
    }
};

}  // namespace detray
