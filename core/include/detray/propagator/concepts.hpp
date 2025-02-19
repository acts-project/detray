/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Include all core actors
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/tuple_helpers.hpp"

// System include(s)
#include <concepts>
#include <optional>

namespace detray::concepts {

/// Concept for a simple actor
template <typename A>
concept actor = std::derived_from<A, detray::actor> &&requires(const A a) {
    typename A::state;
};

/// Concept for an actor including observing actors
template <typename A>
concept composite_actor =
    actor<A> &&A::is_comp_actor::value &&requires(const A ca) {
    typename A::observer_states;
    typename A::observer_state_refs;
};

/// Concept for the actor chain that is being run in the propagator
template <typename A>
concept actor_chain = requires(const A c, typename A::state_tuple s) {
    typename A::actor_tuple;
    typename A::state_tuple;
    typename A::state_ref_tuple;

    { c.actors() }
    ->std::same_as<const typename A::actor_tuple &>;

    { c.make_default_actor_states() }
    ->detray::concepts::any_of<typename A::state_tuple, std::nullopt_t>;

    { c.setup_actor_states(s) }
    ->std::same_as<typename A::state_ref_tuple>;
};

/// Check if a state type belongs to an actor or actor chain
template <typename S, typename T>
concept is_state_of = (actor_chain<T> &&
                       (detail::is_permutation_v<std::remove_cvref_t<S>,
                                                 typename T::state_ref_tuple> ||
                        detail::is_permutation_v<std::remove_cvref_t<S>,
                                                 typename T::state_tuple>)) ||

                      (actor<T> &&
                       std::same_as<std::remove_cvref_t<S>, typename T::state>);

}  // namespace detray::concepts
