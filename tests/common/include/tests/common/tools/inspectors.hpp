/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include "detray/definitions/detail/accessor.hpp"

namespace detray {

/** An inspector that aggregates a number of different inspectors.*/
template <typename... Inspectors>
struct aggregate_inspector {

    using inspector_tuple_t = std::tuple<Inspectors...>;
    inspector_tuple_t _inspectors{};

    template <unsigned int current_id = 0, typename state_type>
    auto operator()(state_type &state, const char *message) {
        // Call inspector
        std::get<current_id>(_inspectors)(state, message);

        // Next mask type
        if constexpr (current_id <
                      std::tuple_size<inspector_tuple_t>::value - 1) {
            return operator()<current_id + 1>(state, message);
        }
    }

    template <typename inspector_t>
    decltype(auto) get() {
        return std::get<inspector_t>(_inspectors);
    }
};

namespace navigation {

/** A navigation inspector that relays information about the encountered
 *  objects the way we need them to compare with the ray
 */
template <int navigation_status = 0,
          template <typename...> class vector_t = dvector>
struct object_tracer {

    // record all object id the navigator encounters
    vector_t<intersection> object_trace = {};

    template <typename state_type>
    auto operator()(state_type &state, const char * /*message*/) {
        // Record the candidate of an encountered object
        if (state.status() == navigation_status) {
            object_trace.push_back(std::move(*(state.current())));
        }
    }

    auto operator[](std::size_t i) { return object_trace[i]; }
};

/** A navigation inspector that prints information about the current navigation
 * state. Meant for debugging.
 */
struct print_inspector {

    // Debug output if an error in the trace is discovered
    std::stringstream debug_stream{};

    template <typename state_type>
    auto operator()(const state_type &state, const char *message) {
        std::string msg(message);
        std::string tabs = "\t\t\t\t\t";

        debug_stream << msg << std::endl;

        debug_stream << "Volume" << tabs << state.volume() << std::endl;
        debug_stream << "surface kernel size\t\t" << state.candidates().size()
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
                debug_stream << " -> " << state.next()->index << std::endl;
            }
        }

        switch (state.status()) {
            case e_abort:
                debug_stream << "status" << tabs << "abort" << std::endl;
                break;
            case e_exit:
                debug_stream << "status" << tabs << "exit" << std::endl;
                break;
            case e_unknown:
                debug_stream << "status" << tabs << "unknowm" << std::endl;
                break;
            case e_towards_object:
                debug_stream << "status" << tabs << "towards_surface"
                             << std::endl;
                break;
            case e_on_target:
                debug_stream << "status" << tabs << "on_surface" << std::endl;
                break;
                /*case e_towards_object:
                    debug_stream << "status" << tabs << "towards_portal" <<
                std::endl; break; case e_on_target: debug_stream << "status" <<
                tabs << "on_portal" << std::endl; break;*/
        };
        debug_stream << "current object\t\t" << state.current_object()
                     << std::endl;
        debug_stream << "distance to next\t";
        if (std::abs(state()) < state.tolerance()) {
            debug_stream << "on obj (within tol)" << std::endl;
        } else {
            debug_stream << state() << std::endl;
        }
        switch (state.trust_level()) {
            case e_no_trust:
                debug_stream << "trust" << tabs << "no_trust" << std::endl;
                break;
            case e_fair_trust:
                debug_stream << "trust" << tabs << "fair_trust" << std::endl;
                break;
            case e_high_trust:
                debug_stream << "trust" << tabs << "high_trust" << std::endl;
                break;
            case e_full_trust:
                debug_stream << "trust" << tabs << "full_trust" << std::endl;
                break;
        };
        debug_stream << std::endl;
    }

    std::string to_string() { return debug_stream.str(); }
};

}  // namespace navigation

namespace propagation {

struct print_inspector {

    std::stringstream stream{};

    template <typename navigator_state_t, typename stepper_state_t>
    DETRAY_HOST_DEVICE void operator()(const navigator_state_t &navigation,
                                       const stepper_state_t &stepping) {

        stream << std::left << std::setw(30);
        switch (static_cast<int>(navigation.status())) {
            case navigation::e_abort:
                stream << "status: abort";
                break;
            case navigation::e_exit:
                stream << "status: exit";
                break;
            case navigation::e_unknown:
                stream << "status: unknowm";
                break;
            case navigation::e_towards_object:
                stream << "status: towards_surface";
                break;
            case navigation::e_on_target:
                stream << "status: on_surface";
                break;
        };

        if (navigation.volume() == dindex_invalid) {
            stream << "volume: " << std::setw(10) << "invalid";
        } else {
            stream << "volume: " << std::setw(10) << navigation.volume();
        }

        if (navigation.current_object() == dindex_invalid) {
            stream << "surface: " << std::setw(14) << "invalid";
        } else {
            stream << "surface: " << std::setw(14)
                   << navigation.current_object();
        }

        stream << "step_size: " << std::setw(10) << stepping._step_size;
    }

    std::string to_string() { return stream.str(); }
};

}  // namespace propagation

namespace stepper {}  // namespace stepper
}  // namespace detray