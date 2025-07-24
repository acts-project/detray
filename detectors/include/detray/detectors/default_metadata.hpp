/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
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
#include "detray/geometry/shapes/unbounded.hpp"
#include "detray/geometry/shapes/unmasked.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/materials/material_map.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/navigation/accelerators/brute_force_finder.hpp"
#include "detray/navigation/accelerators/surface_grid.hpp"

namespace detray {

/// Assembles the detector type. This metatdata contains all available types
template <concepts::algebra algebra_t>
struct default_metadata {

    /// Define the algebra type for the geometry and navigation
    using algebra_type = algebra_t;
    using scalar_t = dscalar<algebra_type>;

    /// Mask-to-(next)-volume link (potentially switchable for SoA)
    using nav_link = std::uint_least16_t;

    /// How to store coordinate transform matrices
    template <template <typename...> class vector_t = dvector>
    using transform_store =
        single_store<dtransform3D<algebra_type>, vector_t, geometry_context>;

    //
    // Surface Primitives
    //

    /// Mask types
    /// @TODO: Need to duplicate for pixel/strip measurement dimensions?
    using rectangle = mask<rectangle2D, algebra_type, nav_link>;
    using trapezoid = mask<trapezoid2D, algebra_type, nav_link>;
    using annulus = mask<annulus2D, algebra_type, nav_link>;
    using cylinder = mask<cylinder2D, algebra_type, nav_link>;
    using cylinder_portal = mask<concentric_cylinder2D, algebra_type, nav_link>;
    using disc = mask<ring2D, algebra_type, nav_link>;
    using straw_tube = mask<line_circular, algebra_type, nav_link>;
    using drift_cell = mask<line_square, algebra_type, nav_link>;
    using single_1 = mask<single3D<1>, algebra_type, nav_link>;
    using single_2 = mask<single3D<2>, algebra_type, nav_link>;
    // TODO: Can single3 be used instead of cylinder portal type or remove it?
    using single_3 = mask<single3D<3>, algebra_type, nav_link>;
    // Debug types, e.g. telescope detector
    using unbounded_rectangle =
        mask<unbounded<rectangle2D>, algebra_type, nav_link>;
    using unbounded_trapezoid =
        mask<unbounded<trapezoid2D>, algebra_type, nav_link>;
    using unbounded_annulus =
        mask<unbounded<annulus2D>, algebra_type, nav_link>;
    using unbounded_cylinder =
        mask<unbounded<cylinder2D>, algebra_type, nav_link>;
    using unbounded_disc = mask<unbounded<ring2D>, algebra_type, nav_link>;
    using unbounded_straw_tube =
        mask<unbounded<line_circular>, algebra_type, nav_link>;
    using unbounded_drift_cell =
        mask<unbounded<line_square>, algebra_type, nav_link>;
    using unmasked_plane = mask<unmasked<>, algebra_type, nav_link>;

    /// Give your mask types a name (needs to be consecutive and has to match
    /// the types position in the mask store!)
    enum class mask_ids : std::uint_least8_t {
        e_rectangle2 = 0u,
        e_portal_rectangle2 = 0u,
        e_trapezoid2 = 1u,
        e_annulus2 = 2u,
        e_cylinder2 = 3u,
        e_portal_cylinder2 = 4u,
        e_ring2 = 5u,
        e_portal_ring2 = 5u,
        e_straw_tube = 6u,
        e_drift_cell = 7u,
        /*e_single1 = 8u,
        e_single2 = 9u,
        e_single3 = 10u,
        e_unbounded_rectangle2 = 11u,
        e_unbounded_trapezoid2 = 12u,
        e_unbounded_annulus2 = 13u,
        e_unbounded_cylinder2 = 14u,
        e_unbounded_disc2 = 15u,
        e_unbounded_line_circular2 = 16u,
        e_unbounded_cell2 = 17u,
        e_unmasked2 = 18u,*/
    };

    /// How to store masks
    template <template <typename...> class vector_t = dvector>
    using mask_store = regular_multi_store<mask_ids, empty_context, dtuple,
                                           vector_t, rectangle, trapezoid,
                                           annulus, cylinder, cylinder_portal,
                                           disc, straw_tube, drift_cell /*,
single_1, single_2, single_3, unbounded_rectangle, unbounded_trapezoid,
unbounded_annulus, unbounded_cylinder, unbounded_disc, unbounded_straw_tube,
unbounded_cell, unmasked_plane*/>;

    //
    // Material Description
    //

    /// Homogeneous material description
    /// @{
    using slab = material_slab<scalar_t>;
    using rod = material_rod<scalar_t>;
    /// @}

    /// Material grid types (default: closed bounds, regular binning)
    /// @{

    // Disc material grid
    template <typename container_t>
    using disc_map_t = material_map<algebra_type, ring2D, container_t>;

    // Cylindrical material grid
    template <typename container_t>
    using cylinder2_map_t = material_map<algebra_type, cylinder2D, container_t>;

    // Concentric cylindrical material grid
    template <typename container_t>
    using concentric_cylinder2_map_t =
        material_map<algebra_type, concentric_cylinder2D, container_t>;

    // Rectangular material grid
    template <typename container_t>
    using rectangular_map_t =
        material_map<algebra_type, rectangle2D, container_t>;

    // Cuboid volume material grid
    template <typename container_t>
    using cuboid_map_t = material_map<algebra_type, cuboid3D, container_t>;

    // Cylindrical volume material grid
    template <typename container_t>
    using cylinder3_map_t = material_map<algebra_type, cylinder3D, container_t>;

    /// @}

    /// Give your material types a name (needs to be consecutive and has to
    /// match the types position in the mask store!)
    enum class material_ids : std::uint_least8_t {
        // Material texture (grid) shapes
        e_concentric_cylinder2_map = 0u,
        e_disc2_map = 1u,
        e_cylinder2_map = 2u,
        e_rectangle2_map = 3u,
        e_trapezoid2_map = 3u,
        e_annulus2_map = 1u,
        e_drift_cell_map = 5u,
        e_straw_tube_map = 5u,
        // Volume material
        e_cuboid3_map = 7,
        e_cylinder3_map = 8u,
        // Homogeneous material
        e_slab = 4u,
        e_rod = 5u,
        e_raw_material = 6u,
        e_none = 9u,
    };

    DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                       material_ids mid) {

        switch (mid) {
            case material_ids::e_concentric_cylinder2_map:
                os << "e_concentric_cylinder2_map";
                break;
            case material_ids::e_disc2_map:
                // e_annulus2_map has same value (1u)
                os << "e_disc2_map/e_annulus2_map";
                break;
            case material_ids::e_cylinder2_map:
                os << "e_cylinder2_map";
                break;
            case material_ids::e_rectangle2_map:
                // e_trapezoid2_map has same value (3u)
                os << "e_rectangle2_map/e_trapezoid2_map";
                break;
            case material_ids::e_slab:
                os << "e_slab";
                break;
            case material_ids::e_drift_cell_map:
                // e_straw_tube_map and e_rod have same value (5u)
                os << "e_drift_cell_map/e_straw_tube_map/e_rod";
                break;
            case material_ids::e_raw_material:
                os << "e_raw_material";
                break;
            case material_ids::e_cuboid3_map:
                os << "e_cuboid3_map";
                break;
            case material_ids::e_cylinder3_map:
                os << "e_cylinder3_map";
                break;
            case material_ids::e_none:
                os << "e_none";
                break;
        }
        return os;
    }

    /// How to store materials
    template <typename container_t = host_container_types>
    using material_store = multi_store<
        material_ids, empty_context, dtuple,
        grid_collection<concentric_cylinder2_map_t<container_t>>,
        grid_collection<disc_map_t<container_t>>,
        grid_collection<cylinder2_map_t<container_t>>,
        grid_collection<rectangular_map_t<container_t>>,
        typename container_t::template vector_type<slab>,
        typename container_t::template vector_type<rod>,
        typename container_t::template vector_type<material<scalar_t>>,
        grid_collection<cuboid_map_t<container_t>>,
        grid_collection<cylinder3_map_t<container_t>>>;

    //
    // Acceleration structures
    //

    /// surface grid types (default boundaries: closed binning)
    /// @TODO: Will we need the types for all grid configurations (binnning,
    /// bin boundaries, serializers)?
    /// @{

    // surface grid definition: bin-content: darray<surface_type, 9>
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

    // 3D cylindrical grid for the entire barrel / endcaps
    template <typename bin_entry_t, typename container_t>
    using cylinder3D_sf_grid =
        surface_grid_t<axes<cylinder3D>, bin_entry_t, container_t>;

    // Irregular binning (hopefully not needed)

    // cylindrical grid for the barrel layers
    template <typename bin_entry_t, typename container_t>
    using irr_cylinder2D_sf_grid =
        surface_grid_t<axes<cylinder2D, axis::bounds::e_closed, axis::irregular,
                            axis::irregular>,
                       bin_entry_t, container_t>;

    // disc grid for the endcap layers
    template <typename bin_entry_t, typename container_t>
    using irr_disc_sf_grid = surface_grid_t<
        axes<ring2D, axis::bounds::e_closed, axis::irregular, axis::irregular>,
        bin_entry_t, container_t>;

    // 3D cylindrical grid for the entire barrel
    template <typename bin_entry_t, typename container_t>
    using irr_cylinder3D_sf_grid =
        surface_grid_t<axes<cylinder3D, axis::bounds::e_closed, axis::irregular,
                            axis::irregular, axis::irregular>,
                       bin_entry_t, container_t>;

    /// @}

    /// Acceleration data structures
    enum class accel_ids : std::uint_least8_t {
        e_brute_force = 0u,     // test all surfaces in a volume (brute force)
        e_cylinder2_grid = 1u,  // e.g. barrel layers
        e_disc_grid = 2u,       // e.g. endcap layers
        e_irr_cylinder2_grid = 3u,
        e_irr_disc_grid = 4u,
        // e_cylinder3_grid = 5u,
        // e_irr_cylinder3_grid = 6u,
        // ... e.g. frustum navigation types
        e_default = e_brute_force,
    };

    DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                       accel_ids aid) {

        switch (aid) {
            case accel_ids::e_brute_force:
                os << "e_brute_force";
                break;
            case accel_ids::e_cylinder2_grid:
                os << "e_cylinder2_grid";
                break;
            case accel_ids::e_disc_grid:
                os << "e_disc_grid";
                break;
            case accel_ids::e_irr_cylinder2_grid:
                os << "e_irr_cylinder2_grid";
                break;
            case accel_ids::e_irr_disc_grid:
                os << "e_irr_disc_grid";
                break;
        }
        return os;
    }

    /// How to link to the entries in the data stores
    using transform_link = typename transform_store<>::link_type;
    using mask_link = typename mask_store<>::range_link;
    using material_link = typename material_store<>::single_link;
    /// Surface type used for sensitives, passives and portals
    using surface_type =
        surface_descriptor<mask_link, material_link, transform_link, nav_link>;

    /// How to store the acceleration data structures
    template <typename container_t = host_container_types>
    using accelerator_store =
        multi_store<accel_ids, empty_context, dtuple,
                    brute_force_collection<surface_type, container_t>,
                    grid_collection<
                        cylinder2D_sf_grid<surface_type, container_t>>,
                    grid_collection<disc_sf_grid<surface_type, container_t>>,
                    grid_collection<
                        irr_cylinder2D_sf_grid<surface_type, container_t>>,
                    grid_collection<
                        irr_disc_sf_grid<surface_type, container_t>> /*,
grid_collection<cylinder3D_sf_grid<surface_type,
container_t>>,
grid_collection<irr_cylinder3D_sf_grid<surface_type,
container_t>>*/>;

    //
    // Volume descriptors
    //

    /// How to index the constituent objects (surfaces) in a volume
    /// If they share the same index value here, they will be added into the
    /// same acceleration data structure (brute force is always at 0)
    enum geo_objects : std::uint_least8_t {
        e_portal = 0u,     // Brute force search
        e_sensitive = 1u,  // Grid accelerated search (can be different types)
        e_passive = 0u,    // Brute force search
        e_size = 2u,     // Every volume holds two acceleration data structures
        e_all = e_size,  // i.e. the brute force method and one grid type
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

    /// How a volume links to the accelration data structures
    /// In this case: One link for portals/passives and one sensitive surfaces
    using object_link_type =
        dmulti_index<dtyped_index<accel_ids, dindex>, geo_objects::e_size>;

    //
    // Volume acceleration structure
    //

    /// Volume search grid
    template <typename container_t = host_container_types>
    using volume_finder =
        grid<algebra_type,
             axes<cylinder3D, axis::bounds::e_open, axis::irregular,
                  axis::regular, axis::irregular>,
             bins::single<dindex>, simple_serializer, container_t>;
};

}  // namespace detray
