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

/// Base class actor implementation
struct actor {
    /// Tag whether this is a composite type
    struct is_comp_actor : public std::false_type {};

    /// Defines the actors state
    struct base_state {};

    // Implicit tag as actor that is inherited by all actors (needed only to
    // compile the conditional in composite_actor)
    struct actor_type {
        using type = actor;
        using state_type = base_state;
    };
};

/// Composition of actors
///
/// The composition represents an actor together with its observers. In
/// addition to running its own implementation, it notifies its observing actors
///
/// @tparam tuple_t the tuple used to unroll observer types.
/// @tparam actor_impl_t the actor the compositions implements itself.
/// @tparam observers a pack of observing actors that get called on the updated
///         actor state of the compositions actor implementation.
template <template <typename...> class tuple_t = dtuple,
          class actor_impl_t = actor, typename... observers>
class composite_actor : public actor_impl_t {

    public:
    /// Tag whether this is a composite type (hides the def in the actor)
    struct is_comp_actor : public std::true_type {};

    /// The composite is an actor in itself. If it is derived from another
    /// composition, it will just implement the same actor as its base class.
    /// I.e. it is impossible to observe another composition's observers.
    using actor_type =
        std::conditional_t<static_cast<bool>(
                               typename actor_impl_t::is_comp_actor()),
                           typename actor_impl_t::actor_type, actor_impl_t>;
    using state_type = typename actor_type::state_type;

    /// Call to the implementation of the actor (the actor possibly being an
    /// observer itself)
    ///
    /// First runs its own implementation, then passes the updated state to its
    /// observers.
    ///
    /// @param states the states of all actors in the chain
    /// @param p_state the state of the propagator (stepper and navigator)
    /// @param subject_state the state of the actor this actor observes. Uses
    ///                      a dummy type if this is not an observing actor.
    template <typename actor_states_t, typename propagator_state_t,
              typename subj_state_t = typename actor::base_state>
    DETRAY_HOST_DEVICE void operator()(
        actor_states_t &states, propagator_state_t &p_state,
        subj_state_t &&subject_state = {}) const {
        // Do your own work ...
        // Two cases: simple actor or observing actor (needs subject state)
        if constexpr (std::is_same_v<subj_state_t,
                                     typename actor::base_state>) {
            static_cast<actor_type const *>(this)->operator()(
                detail::get<typename actor_type::state_type &>(states),
                p_state);
        } else {
            static_cast<actor_type const *>(this)->operator()(
                detail::get<typename actor_type::state_type &>(states), p_state,
                subject_state);
        }

        // Then run the observers on the updated state
        notify(_observers, states,
               detail::get<typename actor_type::state_type &>(states), p_state,
               std::make_index_sequence<sizeof...(observers)>{});
    }

    private:
    /// Notifies the observing actors
    ///
    /// In order to distinguish between actors and composite actors, the
    /// template signature is resolved.
    ///
    /// @param observer one of the observers
    /// @param states the states of all actors in the chain
    /// @param actor_state the state of this compositions actor as the subject
    ///                    to all of its observers
    /// @param p_state the state of the propagator (stepper and navigator)
    template <typename observer_t, typename actor_states_t,
              typename actor_impl_state_t, typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void notify(const observer_t &observer,
                                          actor_states_t &states,
                                          actor_impl_state_t &actor_state,
                                          propagator_state_t &p_state) const {
        // Two cases: observer is a simple actor or a composite actor
        if constexpr (not typename observer_t::is_comp_actor()) {
            observer(detail::get<typename observer_t::state_type &>(states),
                     actor_state, p_state);
        } else {
            observer(states, actor_state, p_state);
        }
    }

    /// Resolve the observer notification.
    ///
    /// Unrolls the observer types and runs the notification for each of them.
    ///
    /// @param observer_list all observers of the actor
    /// @param states the states of all actors in the chain
    /// @param actor_state the state of this compositions actor as the subject
    ///                    to all of its observers
    /// @param p_state the state of the propagator (stepper and navigator)
    template <std::size_t... indices, typename actor_states_t,
              typename actor_impl_state_t, typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void notify(
        const tuple_t<observers...> &observer_list, actor_states_t &states,
        actor_impl_state_t &actor_state, propagator_state_t &p_state,
        std::index_sequence<indices...> /*ids*/) const {

        (notify(detail::get<indices>(observer_list), states, actor_state,
                p_state),
         ...);
    }

    /// Keep the observers (might be composites again)
    tuple_t<observers...> _observers = {};
};

}  // namespace detray
