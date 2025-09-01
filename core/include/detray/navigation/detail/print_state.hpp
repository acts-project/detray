/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/navigation/detail/print_state.hpp"
#include "detray/navigation/navigation_config.hpp"

// System include(s)
#include <iomanip>
#include <sstream>
#include <string>

namespace detray::navigation {

/// Print basic information about the state of a navigation stream @param state
template <typename state_type>
DETRAY_HOST inline std::string print_state(const state_type &state) {

    // Gathers navigation information accross navigator update calls
    std::stringstream debug_stream{};
    // Column width in output
    constexpr int cw{20};

    debug_stream << std::left << std::setw(cw) << "Volume:" << state.volume()
                 << std::endl;

    debug_stream << std::setw(cw) << std::boolalpha
                 << "hearbeat:" << state.is_alive() << std::endl;
    std::cout << std::noboolalpha;

    // Navigation direction
    debug_stream << std::setw(cw) << "direction:";
    switch (state.direction()) {
        using enum direction;
        case e_backward:
            debug_stream << "backward";
            break;
        case e_forward:
            debug_stream << "forward";
            break;
        default:
            break;
    }
    debug_stream << std::endl;

    // Navigation status
    debug_stream << std::setw(cw) << "status:";
    switch (state.status()) {
        using enum status;
        case e_abort:
            debug_stream << "abort";
            break;
        case e_on_target:
            debug_stream << "on_target";
            break;
        case e_unknown:
            debug_stream << "unknowm";
            break;
        case e_towards_object:
            debug_stream << "towards_object";
            break;
        case e_on_module:
            debug_stream << "on_object";
            break;
        case e_on_portal:
            debug_stream << "on_portal";
            break;
        default:
            break;
    }
    debug_stream << std::endl;

    // Navigation trust level
    debug_stream << std::setw(cw) << "trust:";
    switch (state.trust_level()) {
        using enum trust_level;
        case e_no_trust:
            debug_stream << "no_trust";
            break;
        case e_fair:
            debug_stream << "fair_trust";
            break;
        case e_high:
            debug_stream << "high_trust";
            break;
        case e_full:
            debug_stream << "full_trust";
            break;
        default:
            break;
    }
    debug_stream << std::endl;

    // Number of reachable candidates
    debug_stream << std::setw(cw) << "No. reachable:" << state.n_candidates()
                 << std::endl;

    // Current surface
    debug_stream << std::setw(cw) << "current object:";
    if (state.is_on_surface()) {
        // If "exit" is called twice, the state has been cleared
        debug_stream << state.barcode() << std::endl;
    } else if (state.status() == status::e_on_target) {
        debug_stream << "exited" << std::endl;
    } else {
        debug_stream << "undefined" << std::endl;
    }

    // Next surface
    if (!state.candidates().empty()) {
        debug_stream << std::setw(cw) << "next object:";
        if (state.n_candidates() == 0u) {
            debug_stream << "exhausted" << std::endl;
        } else {
            debug_stream << state.next_surface().barcode() << std::endl;
        }
    }

    // Distance to next
    debug_stream << std::setw(cw) << "distance to next:";
    if (!state.is_exhausted() && state.is_on_surface()) {
        debug_stream << "on obj (within tol)" << std::endl;
    } else if (state.is_exhausted()) {
        debug_stream << "no target" << std::endl;
    } else {
        debug_stream << state() << std::endl;
    }

    return debug_stream.str();
}

/// Print candidate and cofiguration information of a navigation state
///
/// @param state the state object of the navigation stream
/// @param cfg the navigation confuration object
/// @param track_pos the current track position
/// @param track_dir the current track direction
template <typename state_type, concepts::point3D point3_t,
          concepts::vector3D vector3_t>
DETRAY_HOST inline std::string print_candidates(const state_type &state,
                                                const navigation::config &cfg,
                                                const point3_t &track_pos,
                                                const vector3_t &track_dir) {

    using detector_t = typename state_type::detector_type;
    using geo_ctx_t = typename detector_t::geometry_context;
    using scalar_t = typename detector_t::scalar_type;

    // Gathers navigation information accross navigator update calls
    std::stringstream debug_stream{};
    // Column width in output
    constexpr int cw{20};

    debug_stream << std::left << std::setw(cw) << "Overstep tol.:"
                 << cfg.overstep_tolerance / detray::unit<scalar_t>::um << " um"
                 << std::endl;

    debug_stream << std::setw(cw) << "Track:"
                 << "pos: [r = " << vector::perp(track_pos)
                 << ", z = " << track_pos[2] << "]," << std::endl;

    debug_stream << std::setw(cw) << " "
                 << "dir: [" << track_dir[0] << ", " << track_dir[1] << ", "
                 << track_dir[2] << "]" << std::endl;

    debug_stream << "Surface candidates: " << std::endl;

    for (const auto &sf_cand : state) {

        debug_stream << std::left << std::setw(6) << "-> " << sf_cand;

        assert(!sf_cand.sf_desc.barcode().is_invalid());

        // Use additional debug information that was gathered on the cand.
        if constexpr (state_type::value_type::contains_pos()) {
            const auto &local = sf_cand.local();
            if (!sf_cand.sf_desc.barcode().is_invalid()) {
                point3_t pos =
                    geometry::surface{state.detector(), sf_cand.sf_desc}
                        .local_to_global(geo_ctx_t{}, local, track_dir);
                debug_stream << " glob: [r = " << vector::perp(pos)
                             << ", z = " << pos[2] << "]" << std::endl;
            } else {
                debug_stream << "Invalid barcode" << std::endl;
            }
        } else {
            debug_stream << std::endl;
        }
    }

    return debug_stream.str();
}

}  // namespace detray::navigation
