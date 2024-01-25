/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/utils/tuple_helpers.hpp"

// System include(s)
#include <iomanip>
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
    auto operator()(state_type &state, const navigation::config &cfg,
                    const char *message) {
        // Call inspector
        std::get<current_id>(_inspectors)(state, cfg, message);

        // Next inspector
        if constexpr (current_id <
                      std::tuple_size<inspector_tuple_t>::value - 1) {
            return operator()<current_id + 1>(state, cfg, message);
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
template <typename candidate_t, template <typename...> class vector_t = dvector,
          status... navigation_status>
struct object_tracer {

    // record all object id the navigator encounters
    vector_t<candidate_t> object_trace = {};

    /// Inspector interface
    template <typename state_type>
    auto operator()(state_type &state, const navigation::config &,
                    const char * /*message*/) {

        // Record the candidate of an encountered object
        if ((is_status(state.status(), navigation_status) or ...)) {
            object_trace.push_back(std::move(*(state.current())));
        }
    }

    /// @returns a specific candidate from the trace
    const candidate_t &operator[](std::size_t i) const {
        return object_trace[i];
    }

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
    auto operator()(const state_type &state, const navigation::config &cfg,
                    const char *message) {
        std::string msg(message);
        std::string tabs = "\t\t\t\t";

        debug_stream << msg << std::endl;

        debug_stream << "Volume" << tabs << state.volume() << std::endl;
        debug_stream << "No. reachable\t\t\t" << state.n_candidates()
                     << std::endl;

        debug_stream << "Surface candidates: " << std::endl;

        using geo_ctx_t = typename state_type::detector_type::geometry_context;
        for (const auto &sf_cand : state.candidates()) {
            const auto &local = sf_cand.local;
            const auto pos =
                surface{*state.detector(), sf_cand.sf_desc}.local_to_global(
                    geo_ctx_t{}, local, {});

            debug_stream << sf_cand;
            debug_stream << ", glob: [r:" << getter::perp(pos)
                         << ", z:" << pos[2] << "]" << std::endl;
        }
        if (not state.candidates().empty()) {
            debug_stream << "=> next: ";
            if (state.is_exhausted()) {
                debug_stream << "exhausted" << std::endl;
            } else {
                debug_stream << " -> " << state.next_surface().barcode()
                             << std::endl;
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

        debug_stream << "current object\t\t\t";
        if (state.is_on_portal() or state.is_on_module() or
            state.status() == status::e_on_target) {
            debug_stream << state.barcode() << std::endl;
        } else {
            debug_stream << "undefined" << std::endl;
        }

        debug_stream << "distance to next\t\t";
        if (math::abs(state()) < cfg.on_surface_tolerance) {
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

namespace stepping {

/// A stepper inspector that prints information about the current stepper
/// state. Meant for debugging.
struct print_inspector {

    /// Gathers stepping information from inside the stepper methods
    std::stringstream debug_stream{};

    /// Inspector interface. Gathers detailed information during stepping
    template <typename state_type>
    void operator()(const state_type &state, const stepping::config &,
                    const char *message) {
        std::string msg(message);
        std::string tabs = "\t\t\t\t";

        debug_stream << msg << std::endl;

        debug_stream << "Step size" << tabs << state.step_size() << std::endl;
        debug_stream << "Path length" << tabs << state.path_length()
                     << std::endl;

        switch (state.direction()) {
            case step::direction::e_forward:
                debug_stream << "direction" << tabs << "forward" << std::endl;
                break;
            case step::direction::e_unknown:
                debug_stream << "direction" << tabs << "unknown" << std::endl;
                break;
            case step::direction::e_backward:
                debug_stream << "direction" << tabs << "backward" << std::endl;
                break;
        };

        auto pos = state().pos();

        debug_stream << "Pos:\t[r = " << math::hypot(pos[0], pos[1])
                     << ", z = " << pos[2] << "]" << std::endl;
        debug_stream << "Tangent:\t"
                     << detail::ray<__plugin::transform3<scalar>>(state())
                     << std::endl;
        debug_stream << std::endl;
    }

    /// Inspector interface. Gathers detailed information during stepping
    template <typename state_type>
    void operator()(const state_type &state, const stepping::config &,
                    const char *message, const std::size_t n_trials,
                    const scalar step_scalor) {
        std::string msg(message);
        std::string tabs = "\t\t\t\t";

        debug_stream << msg << std::endl;

        // Remove trailing newlines
        debug_stream << "Step size" << tabs << state.step_size() << std::endl;
        debug_stream << "no. RK adjustments"
                     << "\t\t" << n_trials << std::endl;
        debug_stream << "Step size scale factor"
                     << "\t\t" << step_scalor << std::endl;

        debug_stream << "Bfield points:" << std::endl;
        const auto &f = state._step_data.b_first;
        debug_stream << "\tfirst:" << tabs << f[0] << ", " << f[1] << ", "
                     << f[2] << std::endl;
        const auto &m = state._step_data.b_middle;
        debug_stream << "\tmiddle:" << tabs << m[0] << ", " << m[1] << ", "
                     << m[2] << std::endl;
        const auto &l = state._step_data.b_last;
        debug_stream << "\tlast:" << tabs << l[0] << ", " << l[1] << ", "
                     << l[2] << std::endl;

        debug_stream << std::endl;
    }

    /// @returns a string representation of the gathered information
    std::string to_string() { return debug_stream.str(); }
};

}  // namespace stepping

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

        if (detail::is_invalid_value(navigation.volume())) {
            printer.stream << "volume: " << std::setw(10) << "invalid";
        } else {
            printer.stream << "volume: " << std::setw(10)
                           << navigation.volume();
        }

        printer.stream << "surface: " << std::setw(14);
        if (navigation.is_on_portal() or navigation.is_on_module()) {
            printer.stream << navigation.barcode();
        } else {
            printer.stream << "undefined";
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
