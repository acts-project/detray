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

/// Wraps the state of an actor
template <typename actor_t>
struct result {
    using state_type = typename actor_t::state;

    state_type &state() { return _state; }
    const state_type &state() const { return _state; }

    state_type _state{};
};

/// Base class actor implementation
///
/// It implements the call and observer handling infrastructure that all
/// actors need to comply with.
///
struct actor {

    /// Defines the actors state
    struct state {};

    // Nothing to be notified
};

/// Composition of actors
///
/// The composite represents an actor together with its observers. It delegates
/// all calls to its observer actors in addition to its own implementation.
///
/// @tparam actor_t the type of actor implemented by the composition.
/// @tparam observers the component actor that are the observers to this actor.
template <template <typename...> class tuple_t = dtuple,
          typename actor_impl_t = actor, typename... observers>
struct composite_actor {

    public:
    using composite_type = composite_actor<tuple_t, actor_impl_t, observers...>;
    /// The composite is an actor in itself
    using actor_type = actor_impl_t;
    /// Aggregate this actors state with its component states
    using actor_result_type = result<actor_type>;
    /// Register the components as observers
    using observer_list_type = tuple_t<observers...>;
    /// List of observer result types
    using observer_results_type = tuple_t<result<observers>...>;
    /// Complete result of the composite
    // using result_type = std::pair<actor_result_type, observer_results_type>;

    /// Contains the observers states
    struct state : public actor_impl_t::state {

        /// Additional data from components
        observer_results_type observer_states{};
    };

    /// Do this actors work and then call the observers on the updated status
    template <typename subj_result_t>
    void operator()(typename result<composite_type>::state_type &comp_state,
                    subj_result_t &subject_result) {
        // Do your own work ...
        _actor(comp_state, subject_result);

        // Then run the observers on the updated state
        notify(comp_state, _observers, comp_state.observer_states,
               std::make_index_sequence<sizeof...(observers)>{});
    }

    /// Do this actors work and then call the observers on the updated status
    void operator()(typename result<composite_type>::state_type &comp_state) {
        // Do your own work ...
        _actor(comp_state);

        // Then run the observers on the updated state
        notify(comp_state, _observers, comp_state.observer_states,
               std::make_index_sequence<sizeof...(observers)>{});
    }

    private:
    /// Call all observers on this actors result
    template <typename state_t, typename observer_t, typename observer_result_t>
    constexpr inline void notify(state_t &comp_state, observer_t &observer,
                                 observer_result_t &observer_result) {
        observer(observer_result.state(), comp_state);
    }

    /// Call all observers on this actors result
    template <typename state_t, std::size_t... indices>
    constexpr inline void notify(
        state_t &comp_state, tuple_t<observers...> &observer_list,
        tuple_t<result<observers>...> &observer_results,
        std::index_sequence<indices...> /*ids*/) {
        (notify(comp_state, detail::get<indices>(observer_list),
                detail::get<indices>(observer_results)),
         ...);
    }

    /// Keep the observers (might be composites again)
    actor_type _actor = {};
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
    /// Aggregated type of all actors results
    using actor_results_type = tuple_t<result<actors_t>...>;

    /// State of the actor chain
    struct state {
        // Aggregates all result types from the registered actors
        using state_type = actor_results_type;

        state_type _data{};
    };

    /// Do this actors work and then call the observers on the updated status
    void operator()(state &results) {

        // Run all actors and their observers
        run(_actors, results._data,
            std::make_index_sequence<sizeof...(actors_t)>{});
    }

    private:
    /// Call all observers on this actors result
    template <typename actor_t>
    constexpr inline void run(actor_t &actor, result<actor_t> &actor_result) {
        actor(actor_result.state());
    }

    /// Call all observers on this actors result
    template <std::size_t... indices>
    constexpr inline void run(actor_list_type &actors,
                              actor_results_type &actor_results,
                              std::index_sequence<indices...> /*ids*/) {
        (run(detail::get<indices>(actors), detail::get<indices>(actor_results)),
         ...);
    }

    private:
    actor_list_type _actors = {};

    /// Get the index for a type. Needs to be unrolled in case of thrust tuple.
    /*template <typename actor_t>
    DETRAY_HOST_DEVICE static constexpr unsigned int get_actor() {
        return unroll_ids<std::remove_reference_t<object_t>,
                          registered_types...>();
    }

    /// Checks whether a given types is known in the registry.
    template <typename object_t>
    DETRAY_HOST_DEVICE static constexpr bool is_registered() {
        return not(get_id<object_t>() == e_unknown);
    }

    /// Extract an index and check it.
    template <typename object_t>
    struct get_index {
        static constexpr ID value = get_id<object_t>();
        constexpr bool operator()() noexcept { return is_valid(value); }
    };

    /// Convert index to ID and do some (limited) checking.
    ///
    /// @tparam ref_idx matches to index arg to perform static checks
    /// @param index argument to be converted to valid id type
    ///
    /// @return the matching ID type.
    template <std::size_t ref_idx = 0>
    DETRAY_HOST_DEVICE static constexpr ID to_id(const std::size_t index) {
        if (ref_idx == index) {
            // Produce a more helpful error than the usual tuple index error
            static_assert(
                is_valid(ref_idx),
                "Index out of range: Please make sure that indices and type "
                "enums match the number of types in container.");
            return static_cast<ID>(index);
        }
        if constexpr (ref_idx < sizeof...(registered_types) - 1) {
            return to_id<ref_idx + 1>(index);
        }
        // This produces a compiler error when used in type unrolling code
        return static_cast<ID>(sizeof...(registered_types));
    }

    /// Return a type for an index. If the index cannot be mapped, there will be
    /// a compiler error.
    template <ID type_id, template <typename...> class tuple_t = dtuple>
    struct get_type {
        using type = std::remove_reference_t<decltype(
            std::get<type_id>(tuple_t<registered_types...>{}))>;
    };

    private:
    /// dummy type
    struct empty_type {};

    /// Gets the position of a type in a parameter pack, without using tuples.
    template <typename object_t, typename first_t = empty_type,
              typename... remaining_types>
    DETRAY_HOST_DEVICE static constexpr ID unroll_ids() {
        if constexpr (not std::is_same_v<first_t, empty_type> and
                      not std::is_same_v<object_t, first_t>) {
            return unroll_ids<object_t, remaining_types...>();
        }
        if constexpr (std::is_same_v<object_t, first_t>) {
            return n_types - sizeof...(remaining_types) - 1;
        }
        return e_unknown;
    }*/
};

}  // namespace detray