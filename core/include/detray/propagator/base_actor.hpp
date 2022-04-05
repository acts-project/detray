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
///
/// It implements the call and typedef infrastructure that all actors need to
/// comply with. Every actor has an id that defines the position of its state
/// objects in the actor state tuple that is accessible through actor_chain
///
/// @tparam ID the actors id needed to link it to its state
template <std::size_t ID>
class actor {

    public:
    /// The type of the actor is referenced when the actor type of a composition
    /// is resolved
    using actor_type = actor<ID>;

    /// Defines the actors state
    struct state {

        constexpr static std::size_t get_id() { return ID; }
    };

    constexpr static std::size_t get_id() { return ID; }
};

class comp_actor {};

/// Composition of actors
///
/// The composition represents an actor together with its observers. In
/// addition to running its own implementation, it notifies its observing actors
///
/// @tparam ID the id for the composition is the id of the actor it implements.
/// @tparam tuple_t the tuple used to unroll observer types.
/// @tparam actor_impl_t the actor the compositions implements itself.
/// @tparam observers a pack of observing actors that get called on the updated
///         actor state of the compositions actor implementation.
template <std::size_t ID, template <typename...> class tuple_t = dtuple,
          template <std::size_t> class actor_impl_t = actor,
          typename... observers>
class composite_actor : private comp_actor, public actor_impl_t<ID> {

    public:
    /// The composite is an actor in itself. If it is derived from another
    /// composition, it will just implement the same actor as its base class.
    /// I.e. it is impossible to observe another composition's observers.
    using actor_type = typename actor_impl_t<ID>::actor_type;
    /// Register the components as observers
    using observer_list_type = tuple_t<observers...>;

    /// State type for the actor that is implemented by the composition
    struct state : public actor_type::state {};

    /// Call to the implementation of the actor.
    ///
    /// First runs its own implementation, then passes the updated state to its
    /// observers.
    ///
    /// @param states the states of the all actors in the chain
    /// @param p_state the state of the propagator (stepper and navigator)
    template <typename actor_states_t, typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(actor_states_t &states,
                                       propagator_state_t &p_state) const {
        // Do your own work ...
        static_cast<actor_type const *>(this)->operator()(
            detail::get<ID>(states), p_state);

        // Then run the observers on the updated state
        notify(_observers, states, detail::get<ID>(states), p_state,
               std::make_index_sequence<sizeof...(observers)>{});
    }

    /// Call to the implementation of the actor (the actor being an observer
    /// itself)
    ///
    /// First runs its own implementation, then passes the updated state to its
    /// observers.
    ///
    /// @param states the states of all actors in the chain
    /// @param subject_state the state of the actor this actor observes
    /// @param p_state the state of the propagator (stepper and navigator)
    template <typename actor_states_t, typename subj_state_t,
              typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(actor_states_t &states,
                                       subj_state_t &subject_state,
                                       propagator_state_t &p_state) const {
        // Do your own work ...
        static_cast<actor_type const *>(this)->operator()(
            detail::get<ID>(states), subject_state, p_state);

        // Then run the observers on the updated state
        notify(_observers, states, detail::get<ID>(states), p_state,
               std::make_index_sequence<sizeof...(observers)>{});
    }

    private:
    /// Notifies the observing actors, when this actor is a composite itself
    ///
    /// In order to distinguish between actors and composite actors, the
    /// template signature is resolved.
    ///
    /// @param actor_state the state of this compositions actor as the subject
    ///                    to all of its observers
    /// @param observer one of the observers
    /// @param states the states of all actors in the chain
    /// @param p_state the state of the propagator (stepper and navigator)
    template <typename observer_t, typename actor_states_t,
              typename actor_impl_state_t, typename propagator_state_t,
              std::enable_if_t<std::is_base_of_v<comp_actor, observer_t>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline void notify(const observer_t &observer,
                                          actor_states_t &states,
                                          actor_impl_state_t &actor_state,
                                          propagator_state_t &p_state) const {

        observer(states, actor_state, p_state);
    }

    /// Notifies the observing actors, normal case.
    ///
    /// In order to distinguish between actors and composite actors, the
    /// template signature is resolved.
    ///
    /// @param actor_state the state of this compositions actor as the subject
    ///                    to all of its observers
    /// @param observer one of the observers
    /// @param states the states of all actors in the chain
    /// @param p_state the state of the propagator (stepper and navigator)
    template <typename observer_t, typename actor_states_t,
              typename propagator_state_t,
              std::enable_if_t<not std::is_base_of_v<comp_actor, observer_t>,
                               bool> = true>
    DETRAY_HOST_DEVICE inline void notify(
        const observer_t &observer, actor_states_t &states,
        typename actor_type::state &actor_state,
        propagator_state_t &p_state) const {

        observer(detail::get<observer_t::get_id()>(states), actor_state,
                 p_state);
    }

    /// Resolve the observer notification.
    ///
    /// Unrolls the observer types and runs the notification for each of them.
    ///
    /// @param actor_state the state of this compositions actor as the subject
    ///                    to all of its observers
    /// @param observer_list all observers of the actor
    /// @param states the states of all actors in the chain
    /// @param p_state the state of the propagator (stepper and navigator)
    template <std::size_t... indices, typename actor_states_t,
              typename propagator_state_t>
    DETRAY_HOST_DEVICE inline void notify(
        const observer_list_type &observer_list, actor_states_t &states,
        typename actor_type::state &actor_state, propagator_state_t &p_state,
        std::index_sequence<indices...> /*ids*/) const {

        (notify(detail::get<indices>(observer_list), states, actor_state,
                p_state),
         ...);
    }

    /// Keep the observers (might be composites again)
    observer_list_type _observers = {};
};

}  // namespace detray