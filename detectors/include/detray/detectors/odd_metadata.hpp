/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/multi_store.hpp"
#include "detray/core/detail/single_store.hpp"
#include "detray/definitions/detail/containers.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/concentric_cylinder2D.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/geometry/shapes/ring2D.hpp"
#include "detray/materials/material_map.hpp"
#include "detray/navigation/accelerators/brute_force_finder.hpp"
#include "detray/navigation/accelerators/surface_grid.hpp"

// Linear algebra types
#include "detray/definitions/detail/algebra.hpp"

namespace detray {

//
// ODD Detector
//

/// Defines a detector that is design to contain the ODD geometry
struct odd_metadata {

    /// Define the algebra type for the geometry and navigation
    using algebra_type = ALGEBRA_PLUGIN<detray::scalar>;

    /// Portal link type between volumes
    using nav_link = std::uint_least16_t;

    //
    // Surface Primitives
    //

    /// The mask types for the ODD sensitive surfaces
    using rectangle = mask<rectangle2D, nav_link>;
    using trapezoid = mask<trapezoid2D, nav_link>;
    // Types for portals
    using cylinder_portal = mask<concentric_cylinder2D, nav_link>;
    using disc_portal = mask<ring2D, nav_link>;

    /// Assign the mask types to the mask tuple container entries.
    enum class mask_ids : std::uint_least8_t {
        e_rectangle2 = 0u,
        e_trapezoid2 = 1u,
        e_portal_cylinder2 = 2u,
        e_portal_ring2 = 3u,
    };

    /// This is the mask collections tuple (in the detector called 'mask store')
    /// the @c regular_multi_store is a vecemem-ready tuple of vectors of
    /// the detector masks.
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using mask_store =
        regular_multi_store<mask_ids, empty_context, tuple_t, vector_t,
                            rectangle, trapezoid, cylinder_portal, disc_portal>;

    /// How to store and link transforms. The geometry context allows to resolve
    /// the conditions data for e.g. module alignment
    template <template <typename...> class vector_t = dvector>
    using transform_store =
        single_store<dtransform3D<algebra_type>, vector_t, geometry_context>;

    //
    // Material Description
    //
    // Cylindrical material map
    template <typename container_t>
    using cylinder_map_t =
        material_map<concentric_cylinder2D, scalar, container_t>;

    // Disc material map
    template <typename container_t>
    using disc_map_t = material_map<ring2D, scalar, container_t>;

    /// Similar to the mask store, there is a material store, which
    enum class material_ids : std::uint_least8_t {
        e_disc2_map = 0u,
        e_concentric_cylinder2_map = 1u,
        e_none = 4u,
    };

    /// How to store and link materials. The material does not make use of
    /// conditions data ( @c empty_context )
    template <template <typename...> class tuple_t = dtuple,
              typename container_t = host_container_types>
    using material_store = multi_store < material_ids,
          empty_context, tuple_t, grid_collection<disc_map_t<container_t>>,
          grid_collection<cylinder_map_t<container_t>>;

    /// Surface descriptor type used for sensitives, passives and portals
    using transform_link = typename transform_store<>::link_type;
    using mask_link = typename mask_store<>::single_link;
    using material_link = typename material_store<>::single_link;
    using surface_type =
        surface_descriptor<mask_link, material_link, transform_link, nav_link>;

    /// The acceleration data structures live in another tuple that needs to be
    /// indexed correctly:
    enum class accel_ids : std::uint_least8_t {
        e_brute_force = 0,     // test all surfaces in a volume (brute force)
        e_disc_grid = 1,       // endcap
        e_cylinder2_grid = 2,  // barrel
        e_default = e_brute_force,
    };

    /// The tuple store that hold the acceleration data structures for all
    /// volumes.
    template <template <typename...> class tuple_t = dtuple,
              typename container_t = host_container_types>
    using accelerator_store =
        multi_store<accel_ids, empty_context, tuple_t,
                    brute_force_collection<surface_type, container_t>>;

    /// How to index the constituent objects in a volume
    enum geo_objects : std::uint_least8_t {
        e_portal = 0,
        e_passive = 0,
        e_sensitive = 1,
        e_size = 2,
        e_all = e_size,
    };

    /// How a volume finds its constituent objects in the detector containers
    using object_link_type =
        dmulti_index<dtyped_index<accel_ids, dindex>, geo_objects::e_size>;

    /// Data structure that allows to find the current detector volume from a
    /// given position. Here: Uniform grid with a 3D cylindrical shape
    template <typename container_t = host_container_types>
    using volume_finder =
        grid<axes<cylinder3D, axis::bounds::e_open, axis::irregular,
                  axis::regular, axis::irregular>,
             bins::single<dindex>, simple_serializer, container_t>;
};

}  // namespace detray
