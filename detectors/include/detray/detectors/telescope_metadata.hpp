/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/multi_store.hpp"
#include "detray/core/detail/single_store.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/navigation/accelerators/brute_force_finder.hpp"

namespace detray {

/// Defines a telescope detector type with only rectangle portals and one
/// additional kind of contained module surfaces (@tparam mask_shape_t)
template <concepts::algebra algebra_t, typename mask_shape_t = rectangle2D>
struct telescope_metadata {

    /// Define the algebra type for the geometry and navigation
    using algebra_type = algebra_t;
    using scalar_t = dscalar<algebra_type>;

    /// Mask to (next) volume link: next volume(s)
    using nav_link = std::uint_least16_t;

    /// How to store coordinate transform matrices
    template <template <typename...> class vector_t = dvector>
    using transform_store =
        single_store<dtransform3D<algebra_type>, vector_t, geometry_context>;

    //
    // Surface Primitives
    //

    /// Mask types (these types are needed for the portals, which are always
    /// there, and to resolve the wire surface material, i.e. slab vs. rod)
    using rectangle = mask<rectangle2D, algebra_type, nav_link>;
    using straw_tube = mask<line_circular, algebra_type, nav_link>;
    using drift_cell = mask<line_square, algebra_type, nav_link>;

    /// Rectangles are always needed as portals (but the yhave the same type as
    /// module rectangles). Only one additional mask shape is allowed
    enum class mask_ids : std::uint_least8_t {
        e_rectangle2 = 0u,
        e_portal_rectangle2 = 0u,
        e_annulus2 = 1u,
        e_cylinder2 = 1u,
        e_ring2 = 1u,
        e_trapezoid2 = 1u,
        e_single1 = 1u,
        e_single2 = 1u,
        e_single3 = 1u,
        e_straw_tube = 1u,
        e_drift_cell = 1u,
        e_unbounded_annulus2 = 1u,
        e_unbounded_cell2 = 1u,
        e_unbounded_cylinder2 = 1u,
        e_unbounded_disc2 = 1u,
        e_unbounded_rectangle2 = 1u,
        e_unbounded_trapezoid2 = 1u,
        e_unbounded_line_circular2 = 1u,
        e_unmasked2 = 1u,
    };

    /// How to store masks
    template <template <typename...> class vector_t = dvector>
    using mask_store = std::conditional_t<
        std::is_same_v<mask<mask_shape_t, algebra_type, nav_link>, rectangle>,
        regular_multi_store<mask_ids, empty_context, dtuple, vector_t,
                            rectangle>,
        regular_multi_store<mask_ids, empty_context, dtuple, vector_t,
                            rectangle,
                            mask<mask_shape_t, algebra_type, nav_link>>>;

    //
    // Material Description
    //

    /// Material types
    using rod = material_rod<scalar_t>;
    using slab = material_slab<scalar_t>;

    /// Material type ids
    enum class material_ids : std::uint_least8_t {
        e_slab = 0u,
        e_raw_material = 1u,  //< used for homogeneous volume material
        e_rod = 2u,
        e_none = 3u,
    };

    /// How to store materials
    template <typename container_t = host_container_types>
    using material_store = std::conditional_t<
        std::is_same_v<mask<mask_shape_t, algebra_type, nav_link>, drift_cell> |
            std::is_same_v<mask<mask_shape_t, algebra_type, nav_link>,
                           straw_tube>,
        regular_multi_store<material_ids, empty_context, dtuple,
                            container_t::template vector_type, slab,
                            material<scalar_t>, rod>,
        regular_multi_store<material_ids, empty_context, dtuple,
                            container_t::template vector_type, slab,
                            material<scalar_t>>>;

    //
    // Acceleration structures
    //

    /// Acceleration data structures
    enum class accel_ids {
        e_brute_force = 0u,  // test all surfaces in a volume (brute force)
        e_default = e_brute_force,
    };

    /// How to link to the entries in the data stores
    using transform_link = typename transform_store<>::link_type;
    using mask_link = typename mask_store<>::single_link;
    using material_link = typename material_store<>::single_link;
    /// Surface type used for sensitives, passives and portals
    using surface_type =
        surface_descriptor<mask_link, material_link, transform_link, nav_link>;

    /// How to store the brute force search data structure
    template <typename container_t = host_container_types>
    using accelerator_store =
        multi_store<accel_ids, empty_context, dtuple,
                    brute_force_collection<surface_type, container_t>>;

    //
    // Volume descriptors
    //

    /// No grids/other acceleration data structure, everything is brute forced
    enum geo_objects : std::uint_least8_t {
        e_portal = 0u,
        e_sensitive = 1u,
        e_size = 2u,
        e_all = e_size,
    };

    /// One link for all surfaces (in the brute force method)
    using object_link_type =
        dmulti_index<dtyped_index<accel_ids, dindex>, geo_objects::e_size>;

    //
    // Volume acceleration structure
    //

    /// Volume search (only one volume exists)
    template <typename container_t = host_container_types>
    using volume_finder = brute_force_collection<dindex, container_t>;
};

}  // namespace detray
