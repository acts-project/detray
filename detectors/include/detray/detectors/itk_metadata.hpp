/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
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
#include "detray/geometry/shapes/annulus2D.hpp"
#include "detray/geometry/shapes/concentric_cylinder2D.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/geometry/shapes/ring2D.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/materials/material_map.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/navigation/accelerators/brute_force_finder.hpp"
#include "detray/navigation/accelerators/surface_grid.hpp"

// Linear algebra types
#include "detray/definitions/algebra.hpp"

namespace detray {

//
// Detector
//

/// Defines a detector that contains squares, trapezoids and a bounding portal
/// box.
template <concepts::algebra algebra_t>
struct itk_metadata {

    /// Define the algebra type for the geometry and navigation
    using algebra_type = algebra_t;

    /// Portal link type between volumes
    using nav_link = std::uint_least16_t;

    /// How to store and link transforms. The geometry context allows to resolve
    /// the conditions data for e.g. module alignment
    template <template <typename...> class vector_t = dvector>
    using transform_store =
        single_store<dtransform3D<algebra_type>, vector_t, geometry_context>;

    //
    // Surface Primitives
    //

    /// The mask types for the detector sensitive surfaces
    using annulus = mask<annulus2D, algebra_type, nav_link>;
    using rectangle = mask<rectangle2D, algebra_type, nav_link>;

    // Types for portals
    using cylinder_portal = mask<concentric_cylinder2D, algebra_type, nav_link>;
    using disc_portal = mask<ring2D, algebra_type, nav_link>;

    /// Assign the mask types to the mask tuple container entries. It may be a
    /// good idea to have the most common types in the first tuple entries, in
    /// order to minimize the depth of the 'unrolling' before a mask is found
    /// in the tuple
    enum class mask_ids : std::uint_least8_t {
        e_rectangle2 = 0u,
        e_annulus2 = 1u,
        e_portal_cylinder2 = 2u,
        e_portal_ring2 = 3u,
    };

    DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                       mask_ids mid) {

        switch (mid) {
            case mask_ids::e_rectangle2:
                os << "e_rectangle2";
                break;
            case mask_ids::e_annulus2:
                os << "e_annulus2";
                break;
            case mask_ids::e_portal_cylinder2:
                os << "e_portal_cylinder2";
                break;
            case mask_ids::e_portal_ring2:
                os << "e_portal_ring2";
                break;
        }
        return os;
    }

    /// This is the mask collections tuple (in the detector called 'mask store')
    /// the @c regular_multi_store is a vecemem-ready tuple of vectors of
    /// the detector masks.
    template <template <typename...> class vector_t = dvector>
    using mask_store =
        regular_multi_store<mask_ids, empty_context, dtuple, vector_t,
                            rectangle, annulus, cylinder_portal, disc_portal>;

    //
    // Material Description
    //

    /// The material types to be mapped onto the surfaces: Here homogeneous
    /// material
    using slab = material_slab<scalar_t>;

    // Cylindrical material map
    template <typename container_t>
    using cylinder_map_t =
        material_map<algebra_type, concentric_cylinder2D, container_t>;

    // Disc material map
    template <typename container_t>
    using disc_map_t = material_map<algebra_type, ring2D, container_t>;

    /// Similar to the mask store, there is a material store, which
    enum class material_ids : std::uint_least8_t {
        e_concentric_cylinder2_map = 0u,
        e_disc2_map = 1u,
        e_none = 2u,
    };

    DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                       material_ids mid) {

        switch (mid) {
            case material_ids::e_concentric_cylinder2_map:
                os << "e_concentric_cylinder2_map";
                break;
            case material_ids::e_disc2_map:
                os << "e_disc2_map";
                break;
            case material_ids::e_none:
                os << "e_none";
                break;
        }
        return os;
    }

    /// How to store and link materials. The material does not make use of
    /// conditions data ( @c empty_context )
    template <typename container_t = host_container_types>
    using material_store =
        multi_store<material_ids, empty_context, dtuple,
                    grid_collection<cylinder_map_t<container_t>>,
                    grid_collection<disc_map_t<container_t>>>;

    //
    // Acceleration structures
    //

    // surface grid definition: dynamic bin size
    template <typename axes_t, typename bin_entry_t, typename container_t>
    using surface_grid_t =
        grid<algebra_type, axes_t, bins::dynamic_array<bin_entry_t>,
             simple_serializer, container_t, false>;

    // 2D cylindrical grid for the barrel layers
    template <typename bin_entry_t, typename container_t>
    using cylinder2D_sf_grid =
        surface_grid_t<axes<concentric_cylinder2D>, bin_entry_t, container_t>;

    // disc grid for the endcap layers
    template <typename bin_entry_t, typename container_t>
    using disc_sf_grid = surface_grid_t<axes<ring2D>, bin_entry_t, container_t>;

    /// The acceleration data structures live in another tuple that needs to be
    /// indexed correctly:
    enum class accel_ids : std::uint_least8_t {
        e_brute_force = 0u,     // test all surfaces in a volume (brute force)
        e_cylinder2_grid = 1u,  // barrel
        e_disc_grid = 2u,       // endcap
        e_default = e_brute_force,
    };

    DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                       accel_ids aid) {

        switch (aid) {
            case accel_ids::e_brute_force:
                // e_default has same value (0u)
                os << "e_brute_force/e_default";
                break;
            case accel_ids::e_cylinder2_grid:
                os << "e_cylinder2_grid";
                break;
            case accel_ids::e_disc_grid:
                os << "e_disc_grid";
                break;
        }
        return os;
    }

    /// Surface descriptor type used for sensitives, passives and portals
    /// It holds the indices to the surface data in the detector data stores
    /// that were defined above
    using transform_link = typename transform_store<>::link_type;
    using mask_link = typename mask_store<>::range_link;
    using material_link = typename material_store<>::single_link;
    /// Surface type used for sensitives, passives and portals
    using surface_type =
        surface_descriptor<mask_link, material_link, transform_link, nav_link>;

    /// The tuple store that hold the acceleration data structures for all
    /// volumes. Every collection of accelerationdata structures defines its
    /// own container and view type. Does not make use of conditions data
    /// ( @c empty_context )
    template <typename container_t = host_container_types>
    using accelerator_store = multi_store<
        accel_ids, empty_context, dtuple,
        brute_force_collection<surface_type, container_t>,
        grid_collection<cylinder2D_sf_grid<surface_type, container_t>>,
        grid_collection<disc_sf_grid<surface_type, container_t>>>;

    //
    // Volume descriptors
    //

    /// How to index the constituent objects in a volume
    /// If they share the same index value here, they will be added into the
    /// same acceleration data structure in every respective volume
    enum geo_objects : std::uint_least8_t {
        e_portal = 0u,
        e_passive = 0u,
        e_sensitive = 1u,
        e_size = 2u,
        e_all = e_size,
    };

    DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                       geo_objects gobj) {

        switch (gobj) {
            case geo_objects::e_portal:
                // e_passive has same value (0u)
                os << "e_portal/e_passive";
                break;
            case geo_objects::e_sensitive:
                os << "e_sensitive";
                break;
            case geo_objects::e_size:
                // e_all has same value (2u)
                os << "e_size/e_all";
                break;
        }
        return os;
    }

    /// How a volume finds its constituent objects in the detector containers
    /// In this case: One range for sensitive/passive surfaces, one for portals
    using object_link_type =
        dmulti_index<dtyped_index<accel_ids, dindex>, geo_objects::e_size>;

    //
    // Volume acceleration structure
    //

    /// Data structure that allows to find the current detector volume from a
    /// given position. Here: Uniform grid with a 3D cylindrical shape
    template <typename container_t = host_container_types>
    using volume_finder =
        grid<algebra_type,
             axes<cylinder3D, axis::bounds::e_open, axis::irregular,
                  axis::regular, axis::irregular>,
             bins::single<dindex>, simple_serializer, container_t>;
};

}  // namespace detray
