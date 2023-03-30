/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/tuple_helpers.hpp"

// System include(s)
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

namespace detray {

/// An inspector that aggregates a number of different inspectors.
template <typename... Inspectors>
struct aggregate_inspector {

    using inspector_tuple_t = std::tuple<Inspectors...>;
    inspector_tuple_t _inspectors{};

    /// Inspector interface
    template <unsigned int current_id = 0, typename state_type>
    auto operator()(state_type &state, const char *message) {
        // Call inspector
        std::get<current_id>(_inspectors)(state, message);

        // Next inspector
        if constexpr (current_id <
                      std::tuple_size<inspector_tuple_t>::value - 1) {
            return operator()<current_id + 1>(state, message);
        }
    }

    /// @returns a specific inspector
    template <typename inspector_t>
    decltype(auto) get() {
        return std::get<inspector_t>(_inspectors);
    }
};

namespace navigation {

/// A navigation inspector that relays information about the encountered
/// objects whenever the navigator reaches one or more status flags
template <template <typename...> class vector_t = dvector,
          status... navigation_status>
struct object_tracer {

    // record all object id the navigator encounters
    vector_t<line_plane_intersection> object_trace = {};

    /// Inspector interface
    template <typename state_type>
    auto operator()(state_type &state, const char * /*message*/) {

        // Record the candidate of an encountered object
        if ((is_status(state.status(), navigation_status) or ...)) {
            object_trace.push_back(std::move(*(state.current())));
        }
    }

    /// @returns a specific candidate from the trace
    auto operator[](std::size_t i) { return object_trace[i]; }

    /// Compares a navigation status with the tracers references
    bool is_status(const status &nav_stat, const status &ref_stat) {
        return (nav_stat == ref_stat);
    }
};

/// A navigation inspector that prints information about the current navigation
/// state. Meant for debugging.
struct print_inspector {

    /// Gathers navigation information accross navigator update calls
    std::stringstream debug_stream{};

    /// Inspector interface. Gathers detailed information during navigation
    template <typename state_type>
    auto operator()(const state_type &state, const char *message) {
        std::string msg(message);
        std::string tabs = "\t\t\t\t";

        debug_stream << msg << std::endl;

        debug_stream << "Volume" << tabs << state.volume() << std::endl;
        debug_stream << "No. reachable\t\t" << state.n_candidates()
                     << std::endl;

        debug_stream << "Surface candidates: " << std::endl;
        for (const auto &sf_cand : state.candidates()) {
            debug_stream << sf_cand.to_string();
        }
        if (not state.candidates().empty()) {
            debug_stream << "=> next: ";
            if (state.is_exhausted()) {
                debug_stream << "exhausted" << std::endl;
            } else {
                debug_stream << " -> " << state.next_object() << std::endl;
            }
        }

        switch (state.status()) {
            case status::e_abort:
                debug_stream << "status" << tabs << "abort" << std::endl;
                break;
            case status::e_on_target:
                debug_stream << "status" << tabs << "e_on_target" << std::endl;
                break;
            case status::e_unknown:
                debug_stream << "status" << tabs << "unknowm" << std::endl;
                break;
            case status::e_towards_object:
                debug_stream << "status" << tabs << "towards_surface"
                             << std::endl;
                break;
            case status::e_on_module:
                debug_stream << "status" << tabs << "on_module" << std::endl;
                break;
            case status::e_on_portal:
                debug_stream << "status" << tabs << "on_portal" << std::endl;
                break;
        };
        debug_stream << "current object\t\t\t" << state.current_object()
                     << std::endl;
        debug_stream << "distance to next\t\t";
        if (std::abs(state()) < state.tolerance()) {
            debug_stream << "on obj (within tol)" << std::endl;
        } else {
            debug_stream << state() << std::endl;
        }
        switch (state.trust_level()) {
            case trust_level::e_no_trust:
                debug_stream << "trust" << tabs << "no_trust" << std::endl;
                break;
            case trust_level::e_fair:
                debug_stream << "trust" << tabs << "fair_trust" << std::endl;
                break;
            case trust_level::e_high:
                debug_stream << "trust" << tabs << "high_trust" << std::endl;
                break;
            case trust_level::e_full:
                debug_stream << "trust" << tabs << "full_trust" << std::endl;
                break;
        };
        debug_stream << std::endl;
    }

    /// @returns a string representation of the gathered information
    std::string to_string() { return debug_stream.str(); }
};

}  // namespace navigation

namespace propagation {

/// Print inspector that runs as actor in the propagation
struct print_inspector : actor {

    struct state {
        std::stringstream stream{};

        std::string to_string() const { return stream.str(); }
    };

    template <typename propagation_state_t>
    DETRAY_HOST_DEVICE void operator()(
        state &printer, const propagation_state_t &prop_state) const {
        const auto &navigation = prop_state._navigation;
        const auto &stepping = prop_state._stepping;

        printer.stream << std::left << std::setw(30);
        switch (navigation.status()) {
            case navigation::status::e_abort:
                printer.stream << "status: abort";
                break;
            case navigation::status::e_on_target:
                printer.stream << "status: e_on_target";
                break;
            case navigation::status::e_unknown:
                printer.stream << "status: unknowm";
                break;
            case navigation::status::e_towards_object:
                printer.stream << "status: towards_surface";
                break;
            case navigation::status::e_on_module:
                printer.stream << "status: on_module";
                break;
            case navigation::status::e_on_portal:
                printer.stream << "status: on_portal";
                break;
        };

        if (navigation.volume() == dindex_invalid) {
            printer.stream << "volume: " << std::setw(10) << "invalid";
        } else {
            printer.stream << "volume: " << std::setw(10)
                           << navigation.volume();
        }

        if (navigation.current_object().is_invalid()) {
            printer.stream << "surface: " << std::setw(14) << "invalid";
        } else {
            printer.stream << "surface: " << std::setw(14)
                           << navigation.current_object();
        }

        printer.stream << "step_size: " << std::setw(10) << stepping._step_size
                       << std::endl;

        printer.stream
            << std::setw(10)
            << detail::ray<
                   typename propagation_state_t::detector_type::transform3>(
                   stepping())
            << std::endl;
    }
};

}  // namespace propagation

}  // namespace detray