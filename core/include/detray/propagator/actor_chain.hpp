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
        static constexpr std::size_t _id{ID};
    };

    /// Call to the implementation of the actor.
    ///
    /// @param actor_state the state of the actor
    /// @param p_state the state of the propagator (stepper and navigator)
    template <typename propagator_state_t>
    void operator()(actor<ID>::state & /*actor_state*/,
                    propagator_state_t & /*p_state*/) {}

    /// Call to the implementation of the actor.
    ///
    /// @param actor_state the state of the actor
    /// @param subject_state the state of the actor this actor is an observer
    ///                      of (if any)
    /// @param p_state the state of the propagator (stepper and navigator)
    template <typename subj_state_t, typename propagator_state_t>
    void operator()(actor<ID>::state & /*actor_state*/,
                    subj_state_t & /*subject_state*/,
                    propagator_state_t & /*p_state*/) {}
};

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
class composite_actor : public actor_impl_t<ID> {

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
    void operator()(actor_states_t &states, propagator_state_t &p_state) {
        // Do your own work ...
        static_cast<actor_type *>(this)->operator()(detail::get<ID>(states),
                                                    p_state);

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
    void operator()(actor_states_t &states, subj_state_t &subject_state,
                    propagator_state_t &p_state) {
        // Do your own work ...
        static_cast<actor_type *>(this)->operator()(detail::get<ID>(states),
                                                    subject_state, p_state);

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
    template <std::size_t obs_ID, template <typename...> class comp_tuple_t,
              template <std::size_t> class comp_actor_impl_t,
              typename... comp_observers,
              template <std::size_t, template <typename...> class,
                        template <std::size_t> class, typename...>
              class observer_t,
              typename actor_states_t, typename propagator_state_t>
    constexpr inline void notify(
        observer_t<obs_ID, comp_tuple_t, comp_actor_impl_t, comp_observers...>
            &observer,
        actor_states_t &states, typename actor_type::state &actor_state,
        propagator_state_t &p_state) {

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
    template <std::size_t obs_ID, template <std::size_t> class observer_t,
              typename actor_states_t, typename propagator_state_t>
    constexpr inline void notify(observer_t<obs_ID> &observer,
                                 actor_states_t &states,
                                 typename actor_type::state &actor_state,
                                 propagator_state_t &p_state) {

        observer(detail::get<obs_ID>(states), actor_state, p_state);
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
    constexpr inline void notify(observer_list_type &observer_list,
                                 actor_states_t &states,
                                 typename actor_type::state &actor_state,
                                 propagator_state_t &p_state,
                                 std::index_sequence<indices...> /*ids*/) {

        (notify(detail::get<indices>(observer_list), states, actor_state,
                p_state),
         ...);
    }

    /// Keep the observers (might be composites again)
    observer_list_type _observers = {};
};

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