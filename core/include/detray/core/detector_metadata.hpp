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
#include "detray/definitions/containers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/cylinder_portal_intersector.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/unbounded.hpp"
#include "detray/masks/unmasked.hpp"
#include "detray/materials/material_map.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/surface_finders/accelerator_grid.hpp"
#include "detray/surface_finders/brute_force_finder.hpp"

namespace detray {

/// Assembles the detector type. This metatdata contains all available types
struct default_metadata {

    /// Mask-to-(next)-volume link (potentially switchable for SoA)
    using nav_link = std::uint_least16_t;

    /// Mask types
    /// @TODO: Need to duplicate for pixel/strip measurement dimensions?
    using rectangle = mask<rectangle2D<>, nav_link>;
    using trapezoid = mask<trapezoid2D<>, nav_link>;
    using annulus = mask<annulus2D<>, nav_link>;
    using cylinder = mask<cylinder2D<>, nav_link>;
    using cylinder_portal =
        mask<cylinder2D<false, cylinder_portal_intersector>, nav_link>;
    using disc = mask<ring2D<>, nav_link>;
    using straw_wire = mask<line<false>, nav_link>;
    using cell_wire = mask<line<true>, nav_link>;
    using single_1 = mask<single3D<1>, nav_link>;
    using single_2 = mask<single3D<2>, nav_link>;
    // TODO: Can single3 be used instead of cylinder portal type or remove it?
    using single_3 = mask<single3D<3>, nav_link>;
    // Debug types, e.g. telescope detector
    using unbounded_rectangle = mask<unbounded<rectangle2D<>>, nav_link>;
    using unbounded_trapezoid = mask<unbounded<trapezoid2D<>>, nav_link>;
    using unbounded_annulus = mask<unbounded<annulus2D<>>, nav_link>;
    using unbounded_cylinder = mask<unbounded<cylinder2D<true>>, nav_link>;
    using unbounded_disc = mask<unbounded<ring2D<>>, nav_link>;
    using unbounded_straw = mask<unbounded<line<false>>, nav_link>;
    using unbounded_cell = mask<unbounded<line<true>>, nav_link>;
    using unmasked_plane = mask<unmasked, nav_link>;

    /// Material types
    using slab = material_slab<detray::scalar>;
    using rod = material_rod<detray::scalar>;

    /// Material grid types (default: closed bounds, regular binning)
    /// @{

    // Disc material grid
    template <typename container_t>
    using disc_map_t =
        material_map<ring2D<>::axes<>, detray::scalar, container_t>;

    // Cylindrical material grid
    template <typename container_t>
    using cylinder2_map_t =
        material_map<cylinder2D<>::axes<>, detray::scalar, container_t>;

    // Rectangular material grid
    template <typename container_t>
    using rectangular_map_t =
        material_map<rectangle2D<>::axes<>, detray::scalar, container_t>;

    // Cuboid volume material grid
    template <typename container_t>
    using cuboid_map_t =
        material_map<cuboid3D<>::axes<>, detray::scalar, container_t>;

    // Cylindrical volume material grid
    template <typename container_t>
    using cylinder3_map_t =
        material_map<cylinder3D::axes<>, detray::scalar, container_t>;

    /// @}

    /// surface grid types (default boundaries: closed binning)
    /// @TODO: Will we need the types for all grid configurations (binnning,
    /// bin boundaries, serializers)?
    /// @{

    // surface grid definition: bin-content: std::array<surface_type, 9>
    template <typename grid_shape_t, typename bin_entry_t, typename container_t>
    using surface_grid_t =
        grid<coordinate_axes<grid_shape_t, false, container_t>,
             bins::static_array<bin_entry_t, 50>, simple_serializer>;

    // 2D cylindrical grid for the barrel layers
    template <typename bin_entry_t, typename container_t>
    using cylinder2D_sf_grid =
        surface_grid_t<cylinder2D<>::axes<>, bin_entry_t, container_t>;

    // disc grid for the endcap layers
    template <typename bin_entry_t, typename container_t>
    using disc_sf_grid =
        surface_grid_t<ring2D<>::axes<>, bin_entry_t, container_t>;

    // 3D cylindrical grid for the entire barrel / endcaps
    template <typename bin_entry_t, typename container_t>
    using cylinder3D_sf_grid =
        surface_grid_t<cylinder3D::axes<>, bin_entry_t, container_t>;

    // Irregular binning (hopefully not needed)

    // cylindrical grid for the barrel layers
    template <typename bin_entry_t, typename container_t>
    using irr_cylinder2D_sf_grid =
        surface_grid_t<cylinder2D<>::axes<n_axis::bounds::e_closed,
                                          n_axis::irregular, n_axis::irregular>,
                       bin_entry_t, container_t>;

    // disc grid for the endcap layers
    template <typename bin_entry_t, typename container_t>
    using irr_disc_sf_grid =
        surface_grid_t<ring2D<>::axes<n_axis::bounds::e_closed,
                                      n_axis::irregular, n_axis::irregular>,
                       bin_entry_t, container_t>;

    // 3D cylindrical grid for the entire barrel
    template <typename bin_entry_t, typename container_t>
    using irr_cylinder3D_sf_grid = surface_grid_t<
        cylinder3D::axes<n_axis::bounds::e_closed, n_axis::irregular,
                         n_axis::irregular, n_axis::irregular>,
        bin_entry_t, container_t>;

    /// @}

    /// How to store coordinate transform matrices
    template <template <typename...> class vector_t = dvector>
    using transform_store = single_store<__plugin::transform3<detray::scalar>,
                                         vector_t, geometry_context>;

    /// Give your mask types a name (needs to be consecutive and has to match
    /// the types position in the mask store!)
    enum class mask_ids : std::uint8_t {
        e_rectangle2 = 0,
        e_portal_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_annulus2 = 2,
        e_cylinder2 = 3,
        e_portal_cylinder2 = 4,
        e_ring2 = 5,
        e_portal_ring2 = 5,
        e_straw_wire = 6,
        e_cell_wire = 7,
        /*e_single1 = 8,
        e_single2 = 9,
        e_single3 = 10,
        e_unbounded_rectangle2 = 11,
        e_unbounded_trapezoid2 = 12,
        e_unbounded_annulus2 = 13,
        e_unbounded_cylinder2 = 14,
        e_unbounded_disc2 = 15,
        e_unbounded_straw2 = 16,
        e_unbounded_cell2 = 17,
        e_unmasked2 = 18,*/
    };

    /// How to store masks
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using mask_store = regular_multi_store<mask_ids, empty_context, tuple_t,
                                           vector_t, rectangle, trapezoid,
                                           annulus, cylinder, cylinder_portal,
                                           disc, straw_wire,
                                           cell_wire /*,
single_1, single_2, single_3, unbounded_rectangle, unbounded_trapezoid,
unbounded_annulus, unbounded_cylinder, unbounded_disc, unbounded_straw,
unbounded_cell, unmasked_plane*/>;

    /// Give your material types a name (needs to be consecutive and has to
    /// match the types position in the mask store!)
    enum class material_ids : std::uint8_t {
        // Material texture (grid) shapes
        e_disc2_map = 0u,
        e_cylinder2_map = 1u,
        e_rectangle2_map = 2u,
        e_trapezoid2_map = 2u,
        e_annulus2_map = 0u,
        e_cell_wire_map = 4u,
        e_straw_wire_map = 4u,
        // Volume material
        // e_cuboid3_map = 3u,
        // e_cylinder3_map = 4u,
        // Homogeneous mapetrial
        e_slab = 3u,
        e_rod = 4u,
        e_none = 5u,
    };

    /// How to store materials
    template <template <typename...> class tuple_t = dtuple,
              typename container_t = host_container_types>
    using material_store =
        multi_store<material_ids, empty_context, tuple_t,
                    grid_collection<disc_map_t<container_t>>,
                    grid_collection<cylinder2_map_t<container_t>>,
                    grid_collection<rectangular_map_t<container_t>>,
                    /*grid_collection<cuboid_map_t<container_t>>,
                    grid_collection<cylinder3_map_t<container_t>>,*/
                    typename container_t::template vector_type<slab>,
                    typename container_t::template vector_type<rod>>;

    /// How to link to the entries in the data stores
    using transform_link = typename transform_store<>::link_type;
    using mask_link = typename mask_store<>::single_link;
    using material_link = typename material_store<>::single_link;
    /// Surface type used for sensitives, passives and portals
    using surface_type =
        surface_descriptor<mask_link, material_link, transform_link, nav_link>;

    /// How to index the constituent objects (surfaces) in a volume
    /// If they share the same index value here, they will be added into the
    /// same acceleration data structure (brute force is always at 0)
    enum geo_objects : std::uint8_t {
        e_portal = 0,     // Brute force search
        e_sensitive = 1,  // Grid accelerated search (can be different types)
        e_passive = 0,    // Brute force search
        e_size = 2,       // Every volume holds two acceleration data structures
        e_all = e_size,   // i.e. the brute force method and one grid type
    };

    /// Acceleration data structures
    enum class accel_ids : std::uint8_t {
        e_brute_force = 0,     // test all surfaces in a volume (brute force)
        e_disc_grid = 1,       // e.g. endcap layers
        e_cylinder2_grid = 2,  // e.g. barrel layers
        e_irr_disc_grid = 3,
        e_irr_cylinder2_grid = 4,
        // e_cylinder3_grid = 5,
        // e_irr_cylinder3_grid = 6,
        // ... e.g. frustum navigation types
        e_default = e_brute_force,
    };

    /// How a volume links to the accelration data structures
    /// In this case: One link for portals/passives and one sensitive surfaces
    using object_link_type =
        dmulti_index<dtyped_index<accel_ids, dindex>, geo_objects::e_size>;

    /// How to store the acceleration data structures
    template <template <typename...> class tuple_t = dtuple,
              typename container_t = host_container_types>
    using accelerator_store =
        multi_store<accel_ids, empty_context, tuple_t,
                    brute_force_collection<surface_type, container_t>,
                    grid_collection<disc_sf_grid<surface_type, container_t>>,
                    grid_collection<
                        cylinder2D_sf_grid<surface_type, container_t>>,
                    grid_collection<
                        irr_disc_sf_grid<surface_type, container_t>>,
                    grid_collection<irr_cylinder2D_sf_grid<
                        surface_type, container_t>> /*,
grid_collection<cylinder3D_sf_grid<surface_type,
container_t>>,
grid_collection<irr_cylinder3D_sf_grid<surface_type,
container_t>>*/>;

    /// Volume search grid
    template <typename container_t = host_container_types>
    using volume_finder =
        grid<coordinate_axes<
                 cylinder3D::axes<n_axis::bounds::e_open, n_axis::irregular,
                                  n_axis::regular, n_axis::irregular>,
                 true, container_t>,
             bins::single<dindex>, simple_serializer>;
};

}  // namespace detray
