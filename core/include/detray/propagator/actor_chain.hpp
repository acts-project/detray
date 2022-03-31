/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <type_traits>
#include <utility>

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/qualifiers.hpp"

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

    /// Call all actors in the chain.
    ///
    /// @param states the states of the actors.
    /// @param p_state the propagation state.
    template <typename actor_states_t, typename propagator_state_t>
    void operator()(actor_states_t &states, propagator_state_t &p_state) {

        run(_actors, states, p_state,
            std::make_index_sequence<sizeof...(actors_t)>{});
    }

    private:
    /// Call a composition of actors.
    ///
    /// @param comp_actr the composite actor
    /// @param states states of all actors (only bare actors)
    /// @param p_state the state of the propagator (stepper and navigator)
    template <std::size_t ID, template <typename...> class comp_tuple_t,
              template <std::size_t> class actor_impl_t, typename... observers,
              template <std::size_t, template <typename...> class,
                        template <std::size_t> class, typename...>
              class actor_t,
              typename actor_states_t, typename propagator_state_t>
    constexpr inline void run(
        actor_t<ID, comp_tuple_t, actor_impl_t, observers...> &comp_actr,
        actor_states_t &states, propagator_state_t &p_state) {
        comp_actr(states, p_state);
    }

    /// Call a single actor.
    ///
    /// @param actr the actor
    /// @param states states of all actors (only bare actors)
    /// @param p_state the state of the propagator (stepper and navigator)
    template <std::size_t ID, template <std::size_t> class actor_t,
              typename actor_states_t, typename propagator_state_t>
    constexpr inline void run(actor_t<ID> &actr, actor_states_t &states,
                              propagator_state_t &p_state) {
        actr(detail::get<ID>(states), p_state);
    }

    /// Resolve the actor calls.
    ///
    /// @param actors list of all actors
    /// @param states states of all actors (only bare actors)
    /// @param p_state the state of the propagator (stepper and navigator)
    template <std::size_t... indices, typename actor_states_t,
              typename propagator_state_t>
    constexpr inline void run(actor_list_type &actors, actor_states_t &states,
                              propagator_state_t &p_state,
                              std::index_sequence<indices...> /*ids*/) {
        (run(detail::get<indices>(actors), states, p_state), ...);
    }

    private:
    /// Tuple of actors
    actor_list_type _actors = {};
};

}  // namespace detray