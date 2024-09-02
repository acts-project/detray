/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <type_traits>
#include <utility>

#include "detray/definitions/detail/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/tuple_helpers.hpp"

namespace detray {

/// The interface to the actors and aborters in the propagation.
///
/// It can hold both simple actors, as well as an actor with its observers.
/// The states of the actors need to be passed to the chain in an external tuple
///
/// @tparam tuple_t tuple type used to resolve the actor types
/// @tparam actors_t the types of the actors in the chain.
template <template <typename...> class tuple_t = dtuple, typename... actors_t>
class actor_chain {

    public:
    /// Types of the actors that are registered in the chain
    using actor_list_type = tuple_t<actors_t...>;
    // Type of states tuple that is used in the propagator
    using state = tuple_t<typename actors_t::state &...>;

    /// Call all actors in the chain.
    ///
    /// @param states the states of the actors.
    /// @param p_state the propagation state.
    template <typename actor_states_t, typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(actor_states_t &states,
                                       propagator_state_t &p_state) const {

        run(states, p_state, std::make_index_sequence<sizeof...(actors_t)>{});
    }

    /// @returns the actor list
    DETRAY_HOST_DEVICE const actor_list_type &actors() const { return _actors; }

    private:
    /// Call the actors. Either single actor or composition.
    ///
    /// @param actr the actor (might be a composite actor)
    /// @param states states of all actors (only bare actors)
    /// @param p_state the state of the propagator (stepper and navigator)
    template <typename actor_t, typename actor_states_t,
              typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void run(const actor_t &actr,
                                       actor_states_t &states,
                                       propagator_state_t &p_state) const {
        if constexpr (!typename actor_t::is_comp_actor()) {
            actr(detail::get<typename actor_t::state &>(states), p_state);
        } else {
            actr(states, p_state);
        }
    }

    /// Resolve the actor calls.
    ///
    /// @param states states of all actors (only bare actors)
    /// @param p_state the state of the propagator (stepper and navigator)
    template <typename actor_states_t, typename propagator_state_t,
              std::size_t... indices>
    DETRAY_HOST_DEVICE inline void run(
        actor_states_t &states, propagator_state_t &p_state,
        std::index_sequence<indices...> /*ids*/) const {
        (run(detail::get<indices>(_actors), states, p_state), ...);
    }

    /// Tuple of actors
    actor_list_type _actors = {};
};

/// Empty actor chain (placeholder)
template <>
class actor_chain<> {

    public:
    /// Empty states replaces a real actor states container
    struct state {};

    /// Call to actors does nothing.
    ///
    /// @param states the states of the actors.
    /// @param p_state the propagation state.
    template <typename actor_states_t, typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(actor_states_t & /*states*/,
                                       propagator_state_t & /*p_state*/) const {
        /*Do nothing*/
    }
};

}  // namespace detray
