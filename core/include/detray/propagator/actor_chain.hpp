/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <iostream>
#include <type_traits>
#include <utility>

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/// Base class actor implementation
///
/// It implements the call and observer handling infrastructure that all
/// actors need to comply with.
///
template <std::size_t ID>
struct actor {

    /// Defines the actors state
    struct state {
        static constexpr std::size_t _id{ID};
    };

    /// Do this actors work and then call the observers on the updated status
    template <typename propagator_state_t>
    void operator()(actor<ID>::state &results, propagator_state_t &p_state) {}

    /// Do this actors work and then call the observers on the updated status
    template <typename subj_result_t, typename propagator_state_t>
    void operator()(actor<ID>::state &results, subj_result_t &subject_result,
                    propagator_state_t &p_state) {}
};

/// Composition of actors
///
/// The composite represents an actor together with its observers. It delegates
/// all calls to its observer actors in addition to its own implementation.
///
/// @tparam actor_t the type of actor implemented by the composition.
/// @tparam observers the component actor that are the observers to this actor.
template <std::size_t ID, template <typename...> class tuple_t = dtuple,
          template <std::size_t> class actor_impl_t = actor,
          typename... observers>
struct composite_actor : public actor_impl_t<ID> {

    public:
    /// The composite is an actor in itself
    using actor_type = typename actor_impl_t<ID>::actor_type;
    /// Register the components as observers
    using observer_list_type = tuple_t<observers...>;

    /// Inherited state type
    struct state : public actor_type::state {};

    /// Do this actors work and then call the observers on the updated status
    template <typename actor_results_t, typename propagator_state_t>
    void operator()(actor_results_t &results, propagator_state_t &p_state) {
        // Do your own work ...
        static_cast<actor_type *>(this)->operator()(
            detail::get<typename actor_type::state>(results), p_state);

        // Then run the observers on the updated state
        notify(detail::get<typename actor_type::state>(results), _observers,
               results, p_state,
               std::make_index_sequence<sizeof...(observers)>{});
    }

    /// Do this actors work and then call the observers on the updated status
    template <typename actor_results_t, typename subj_result_t,
              typename propagator_state_t>
    void operator()(actor_results_t &results, subj_result_t &subject_result,
                    propagator_state_t &p_state) {
        // Do your own work ...

        static_cast<actor_type *>(this)->operator()(
            detail::get<typename actor_type::state>(results), subject_result,
            p_state);

        // Then run the observers on the updated state
        notify(detail::get<typename actor_type::state>(results), _observers,
               results, p_state,
               std::make_index_sequence<sizeof...(observers)>{});
    }

    private:
    /// Call all observers on this actors result
    template <std::size_t obs_ID, template <typename...> class comp_tuple_t,
              template <std::size_t> class comp_actor_impl_t,
              typename... comp_observers,
              template <std::size_t, template <typename...> class,
                        template <std::size_t> class, typename...>
              class observer_t,
              typename observer_results_t, typename propagator_state_t>
    constexpr inline void notify(
        typename actor_type::state &actor_state,
        observer_t<obs_ID, comp_tuple_t, comp_actor_impl_t, comp_observers...>
            &observer,
        observer_results_t &observer_results, propagator_state_t &p_state) {

        observer(observer_results, actor_state, p_state);
    }

    /// Call all observers on this actors result
    template <std::size_t obs_ID, template <std::size_t> class observer_t,
              typename observer_results_t, typename propagator_state_t>
    constexpr inline void notify(typename actor_type::state &actor_state,
                                 observer_t<obs_ID> &observer,
                                 observer_results_t &observer_results,
                                 propagator_state_t &p_state) {

        observer(
            detail::get<typename observer_t<obs_ID>::state>(observer_results),
            actor_state, p_state);
    }

    /// Call all observers on this actors result
    template <std::size_t... indices, typename observer_results_t,
              typename propagator_state_t>
    constexpr inline void notify(typename actor_type::state &actor_state,
                                 observer_list_type &observer_list,
                                 observer_results_t &observer_results,
                                 propagator_state_t &p_state,
                                 std::index_sequence<indices...> /*ids*/) {

        (notify(actor_state, detail::get<indices>(observer_list),
                observer_results, p_state),
         ...);
    }

    /// Keep the observers (might be composites again)
    observer_list_type _observers = {};
};

/// The actor chain is the interface to a composition of actors
///
///
template <template <typename...> class tuple_t = dtuple, typename... actors_t>
class actor_chain {
    public:
    /// Start of the actor call
    using actor_list_type = tuple_t<actors_t...>;

    /// Do this actors work and then call the observers on the updated status
    template <typename actor_results_t, typename propagator_state_t>
    void operator()(actor_results_t &results, propagator_state_t &p_state) {

        // Run all actors and their observers
        run(_actors, results, p_state,
            std::make_index_sequence<sizeof...(actors_t)>{});
    }

    private:
    /// Call all observers on this actors result
    template <std::size_t ID, template <typename...> class comp_tuple_t,
              template <std::size_t> class actor_impl_t, typename... observers,
              template <std::size_t, template <typename...> class,
                        template <std::size_t> class, typename...>
              class actor_t,
              typename actor_results_t, typename propagator_state_t>
    constexpr inline void run(
        actor_t<ID, comp_tuple_t, actor_impl_t, observers...> &actr,
        actor_results_t &results, propagator_state_t &p_state) {
        actr(results, p_state);
    }

    /// Call all observers on this actors result
    template <std::size_t ID, template <std::size_t> class actor_t,
              typename actor_results_t, typename propagator_state_t>
    constexpr inline void run(actor_t<ID> &actr, actor_results_t &results,
                              propagator_state_t &p_state) {
        actr(detail::get<typename actor_t<ID>::state>(results), p_state);
    }

    /// Call all observers on this actors result
    template <std::size_t... indices, typename actor_results_t,
              typename propagator_state_t>
    constexpr inline void run(actor_list_type &actors, actor_results_t &results,
                              propagator_state_t &p_state,
                              std::index_sequence<indices...> /*ids*/) {
        (run(detail::get<indices>(actors), results, p_state), ...);
    }

    private:
    actor_list_type _actors = {};
};

}  // namespace detray