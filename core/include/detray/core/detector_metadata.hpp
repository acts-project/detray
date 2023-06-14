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
#include "detray/geometry/surface.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/cylinder_portal_intersector.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/unbounded.hpp"
#include "detray/masks/unmasked.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/surface_finders/accelerator_grid.hpp"
#include "detray/surface_finders/brute_force_finder.hpp"

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/vector.hpp>

namespace detray {

/// Defines all available types
struct full_metadata {

    /// mask-to-(next)-volume link (potentially switchable for SoA)
    using nav_link = std::uint_least16_t;

    /// mask types
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

    /// material types
    using slab = material_slab<detray::scalar>;
    using rod = material_rod<detray::scalar>;

    /// @TODO material grids
    /// @{
    /// ...
    /// @}

    /// surface grid types (default bounadries: closed binning)
    /// @TODO: Will we need the types for all grid configurations (binnning,
    /// bin boundaries, serializers)?
    /// @{

    // surface grid definition: bin-content: std::array<dindex, 9>
    template <typename grid_shape_t, typename bin_entry_t, typename container_t>
    using regular_surface_grid_t =
        grid<coordinate_axes<grid_shape_t, false, container_t>, bin_entry_t,
             simple_serializer, regular_attacher<9>>;

    // 2D cylindrical grid for the barrel layers
    template <typename bin_entry_t, typename container_t>
    using cylinder2D_sf_grid =
        regular_surface_grid_t<cylinder2D<>::axes<>, bin_entry_t, container_t>;

    // 3D cylindrical grid for the entire barrel / endcaps
    template <typename bin_entry_t, typename container_t>
    using cylinder3D_sf_grid =
        regular_surface_grid_t<cylinder3D::axes<>, bin_entry_t, container_t>;

    // disc grid for the endcap layers
    template <typename bin_entry_t, typename container_t>
    using disc_sf_grid =
        regular_surface_grid_t<ring2D<>::axes<>, bin_entry_t, container_t>;

    // Irregular binning (hopefully not needed)

    // cylindrical grid for the barrel layers
    /*template <typename bin_entry_t, typename container_t>
    using irr_cylinder_sf_grid =
        surface_grid_t<cylinder2D<>::axes<n_axis::bounds::e_closed,
    n_axis::irregular, n_axis::irregular>, bin_entry_t, container_t>;

    // disc grid for the endcap layers
    template <typename bin_entry_t, typename container_t>
    using irr_disc_sf_grid =
    surface_grid_t<ring2D<>::axes<n_axis::bounds::e_closed, n_axis::irregular,
    n_axis::irregular>, bin_entry_t, container_t>;

    // 3D cylindrical grid for the entire barrel
    template <typename bin_entry_t, typename container_t>
    using irr_cylinder3D_sf_grid =
        surface_grid_t<cylinder3D<>::axes<n_axis::bounds::e_closed,
    n_axis::irregular, n_axis::irregular, n_axis::irregular>, bin_entry_t,
    container_t>;*/

    /// @}

    // @TODO: Switch to inhomogenous b-field and read homogeneous b-field also
    // from file.
    using bfield_backend_t =
        covfie::backend::constant<covfie::vector::vector_d<scalar, 3>,
                                  covfie::vector::vector_d<scalar, 3>>;

    /// How to index the constituent objects in a volume
    /// If they share the same index value here, they will be added into the
    /// same container range without any sorting guarantees
    enum geo_objects : std::size_t {
        e_portal = 0,
        // TODO: Put into grids
        e_sensitive = 0,  // Sensitive and passive surfaces will be held in the
        e_passive = 0,    // same acceleration data structure.
        e_size = 2,       // Every volume holds two acceleration data structures
        e_all = e_size,   // One for portals and one for all other surfaces.
    };

    /// How a volume finds its constituent objects in the detector containers
    /// In this case: One range for sensitive/passive surfaces, one for portals
    using object_link_type = dmulti_index<dindex_range, geo_objects::e_size>;

    /// How to store and link transforms
    template <template <typename...> class vector_t = dvector>
    using transform_store = single_store<__plugin::transform3<detray::scalar>,
                                         vector_t, geometry_context>;

    /// Give your mask types a name (needs to be consecutive to be matched
    /// to a type!)
    enum class mask_ids {
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

    /// How to store and link masks
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

    /// Give your material types a name (needs to be consecutive to be matched
    /// to a type!)
    /// @TODO: Add the material grid types for every surface shape
    enum class material_ids {
        e_slab = 0,
        e_rod = 1,
        // ... material map types
        e_none = 2,
    };

    /// How to store and link materials
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using material_store = regular_multi_store<material_ids, empty_context,
                                               tuple_t, vector_t, slab, rod>;

    /// Surface type used for sensitives, passives and portals
    using transform_link = typename transform_store<>::link_type;
    using mask_link = typename mask_store<>::single_link;
    using material_link = typename material_store<>::single_link;
    using source_link = dindex;
    using surface_type = surface<mask_link, material_link, transform_link,
                                 nav_link, source_link>;

    /// Surface finders
    enum class sf_finder_ids {
        e_brute_force = 0,  // test all surfaces in a volume (brute force)
        // e_disc_grid = 1,
        // e_cylinder2_grid = 2,
        // e_cylinder3_grid = 3,
        // e_irr_disc_grid = 4,
        // e_irr_cylinder2_grid = 5,
        // e_irr_cylinder3_grid = 6,
        // ... e.g. frustum navigation types
        e_default = e_brute_force,
    };

    /// How to store and link surface grids
    template <template <typename...> class tuple_t = dtuple,
              typename container_t = host_container_types>
    using surface_finder_store =
        multi_store<sf_finder_ids, empty_context, tuple_t,
                    brute_force_collection<
                        surface_type, container_t> /*,
  grid_collection<disc_sf_grid<surface_type,
  container_t>>,
  grid_collection<cylinder2D_sf_grid<surface_type,
  container_t>>,
  grid_collection<cylinder3D_sf_grid<surface_type,
  container_t>>,
  grid_collection<irr_disc_sf_grid<surface_type,
  container_t>,
  grid_collection<irr_cylinder2D_sf_grid<surface_type,
  container_t>>,
  grid_collection<irr_cylinder3D_sf_grid<surface_type,
  container_t>>*/>;

    /// Volume grid
    template <typename container_t = host_container_types>
    using volume_finder =
        grid<coordinate_axes<
                 cylinder3D::axes<n_axis::bounds::e_open, n_axis::irregular,
                                  n_axis::regular, n_axis::irregular>,
                 true, container_t>,
             dindex, simple_serializer, replacer>;
};

}  // namespace detray
