/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/detectors/odd_metadata.hpp"

namespace detray {

/// Defines the data types needed for the toy detector (same as ODD)
template <concepts::algebra algebra_t>
struct toy_metadata {

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

    /// Mask types
    using rectangle = mask<rectangle2D, algebra_type, nav_link>;
    using trapezoid = mask<trapezoid2D, algebra_type, nav_link>;

    // using cylinder = mask<cylinder2D, algebra_type, nav_link>;  // beampipe
    using cylinder_portal = mask<concentric_cylinder2D, algebra_type, nav_link>;
    using disc_portal = mask<ring2D, algebra_type, nav_link>;

    /// Mask type ids
    enum class mask_ids : std::uint_least8_t {
        e_rectangle2 = 0u,
        e_trapezoid2 = 1u,
        e_portal_cylinder2 = 2u,
        e_portal_ring2 = 3u,
        e_cylinder2 = 2u,
    };

    DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                       mask_ids mid) {

        switch (mid) {
            case mask_ids::e_rectangle2:
                os << "e_rectangle2";
                break;
            case mask_ids::e_trapezoid2:
                os << "e_trapezoid2";
                break;
            case mask_ids::e_portal_cylinder2:
                // e_cylinder2 has same value (2u)
                os << "e_portal_cylinder2/e_cylinder2";
                break;
            case mask_ids::e_portal_ring2:
                os << "e_portal_ring2";
                break;
            default:
                os << "invalid";
        }
        return os;
    }

    /// How to store masks
    template <template <typename...> class vector_t = dvector>
    using mask_store =
        regular_multi_store<mask_ids, empty_context, dtuple, vector_t,
                            rectangle, trapezoid, cylinder_portal, disc_portal>;

    //
    // Material Description
    //

    // Homogeneous material description
    using slab = material_slab<scalar_t>;

    // Cylindrical material grid
    template <typename container_t>
    using cylinder_map_t =
        material_map<algebra_type, concentric_cylinder2D, container_t>;

    // Disc material grid
    template <typename container_t>
    using disc_map_t = material_map<algebra_type, ring2D, container_t>;

    /// Material type ids
    enum class material_ids : std::uint_least8_t {
        e_concentric_cylinder2_map = 0u,
        e_disc2_map = 1u,
        e_slab = 2u,
        e_none = 3u,
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
            case material_ids::e_slab:
                os << "e_slab";
                break;
            case material_ids::e_none:
                os << "e_none";
                break;
            default:
                os << "invalid";
        }
        return os;
    }

    /// How to store materials
    template <typename container_t = host_container_types>
    using material_store =
        multi_store<material_ids, empty_context, dtuple,
                    grid_collection<cylinder_map_t<container_t>>,
                    grid_collection<disc_map_t<container_t>>,
                    typename container_t::template vector_type<slab>>;

    //
    // Acceleration structures
    //

    /// Surface grid types (regular, open binning)
    /// @{

    // Surface grid definition: bin-content: darray<sf_descriptor, 1>
    template <typename axes_t, typename bin_entry_t, typename container_t>
    using surface_grid_t =
        spatial_grid<algebra_type, axes_t, bins::static_array<bin_entry_t, 1>,
                     simple_serializer, container_t, false>;

    // cylindrical grid for the barrel layers
    template <typename bin_entry_t, typename container_t>
    using cylinder_sf_grid =
        surface_grid_t<axes<concentric_cylinder2D>, bin_entry_t, container_t>;

    // disc grid for the endcap layers
    template <typename bin_entry_t, typename container_t>
    using disc_sf_grid = surface_grid_t<axes<ring2D>, bin_entry_t, container_t>;

    /// @}

    /// How to link to the entries in the data stores
    using transform_link = typename transform_store<>::single_link;
    using mask_link = typename mask_store<>::range_link;
    using material_link = typename material_store<>::single_link;
    /// Surface type used for sensitives, passives and portals
    using surface_type =
        surface_descriptor<mask_link, material_link, transform_link, nav_link>;

    //
    // Volume descriptors
    //

    /// Portals and passives in the brute froce search, sensitives in the grids
    enum geo_objects : std::uint_least8_t {
        e_portal = 0u,
        e_passive = 0u,
        e_sensitive = 1u,
        e_volume = 2u,
        e_size = 3u,
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
            case geo_objects::e_volume:
                os << "e_volume";
                break;
            case geo_objects::e_size:
                // e_all has same value (2u)
                os << "e_size/e_all";
                break;
            default:
                os << "invalid";
        }
        return os;
    }

    /// Acceleration data structures
    enum class accel_ids : std::uint_least8_t {
        e_brute_force = 0u,     // test all surfaces in a volume (brute force)
        e_cylinder2_grid = 1u,  // barrel
        e_disc_grid = 2u,       // endcap
        e_volume_cylinder3_grid = 3u,
        e_default = e_brute_force,
        e_default_volume_searcher = e_volume_cylinder3_grid,
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
            case accel_ids::e_volume_cylinder3_grid:
                os << "e_volume_cylinder3_grid/e_default_volume_searcher";
                break;
            default:
                os << "invalid";
        }
        return os;
    }

    /// One link for portals/passives and one sensitive surfaces
    using object_link_type =
        dmulti_index<dtyped_index<accel_ids, dindex>, geo_objects::e_size>;

    template <typename container_t = host_container_types>
    using volume_accelerator =
        spatial_grid<algebra_type,
                     axes<cylinder3D, axis::bounds::e_open, axis::irregular,
                          axis::regular, axis::irregular>,
                     bins::single<dindex>, simple_serializer, container_t,
                     false>;

    /// How to store the acceleration data structures
    template <typename container_t = host_container_types>
    using accelerator_store = multi_store<
        accel_ids, empty_context, dtuple,
        brute_force_collection<surface_type, container_t>,
        grid_collection<cylinder_sf_grid<surface_type, container_t>>,
        grid_collection<disc_sf_grid<surface_type, container_t>>,
        grid_collection<volume_accelerator<container_t>>>;
};

}  // namespace detray
