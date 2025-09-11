/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/grid/detail/concepts.hpp"
#include "detray/utils/ranges/ranges.hpp"

// System include(s)
#include <concepts>

namespace detray::concepts {

/// Concept for a detray geometry acceleration structures
template <class A>
concept accelerator = requires(const A acc) {
    typename A::value_type;

    // Iterate through all contained geometry objects
    { acc.all() } -> ranges::range_of<typename A::value_type>;

    // TODO: In order to require the search method, we need to pass a detector
    // which is auto-deduced
};

/// Concept for a collection of accelerator (data) that can be stored in the
/// detector
template <class AC>
concept accelerator_collection =
    viewable<AC> && bufferable<AC> &&
    requires(const AC accel_coll, unsigned int idx) {
        typename AC::size_type;
        requires concepts::accelerator<typename AC::value_type>;

        { accel_coll[idx] } -> concepts::accelerator;
    };

/// Acceleration structure that contains surfaces (surface descriptors)
/// TODO: Add surface descriptor concept to geometry package
template <class A>
concept surface_accelerator =
    concepts::accelerator<A> &&
    std::same_as<typename A::value_type,
                 surface_descriptor<typename A::value_type::mask_link,
                                    typename A::value_type::material_link,
                                    typename A::value_type::transform_link,
                                    typename A::value_type::navigation_link>>;

/// Acceleration structure that contains volumes (volume indices)
template <class A>
concept volume_accelerator =
    concepts::accelerator<A> && std::same_as<typename A::value_type, dindex>;

/// Brute force acceleration structures
/// @{
template <class A>
concept brute_force_surface_searcher =
    !concepts::grid<A> && concepts::volume_accelerator<A>;

template <class A>
concept brute_force_volume_searcher =
    !concepts::grid<A> && concepts::surface_accelerator<A>;
/// @}

/// Grid based acceleration structures
/// @{
template <class A>
concept surface_grid = concepts::grid<A> && concepts::surface_accelerator<A>;

template <class A>
concept volume_grid = concepts::grid<A> && concepts::volume_accelerator<A>;
/// @}

}  // namespace detray::concepts
