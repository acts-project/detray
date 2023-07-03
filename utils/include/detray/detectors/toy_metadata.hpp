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
#include "detray/materials/material_slab.hpp"
#include "detray/surface_finders/accelerator_grid.hpp"
#include "detray/surface_finders/brute_force_finder.hpp"

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/vector.hpp>

namespace detray {

/// Defines the data types needed for the toy detector
template <typename _bfield_backend_t =
              covfie::backend::constant<covfie::vector::vector_d<scalar, 3>,
                                        covfie::vector::vector_d<scalar, 3>>>
struct toy_metadata {

    /// Mask to (next) volume link: next volume(s)
    using nav_link = std::uint_least16_t;

    /// Mask types
    using rectangle = mask<rectangle2D<>, nav_link>;
    using trapezoid = mask<trapezoid2D<>, nav_link>;
    // using cylinder = mask<cylinder2D<>, nav_link>;  // beampipe
    using cylinder_portal =
        mask<cylinder2D<false, cylinder_portal_intersector>, nav_link>;
    using disc_portal = mask<ring2D<>, nav_link>;

    /// Material types
    using slab = material_slab<detray::scalar>;

    /// surface grid types (regular, open binning)
    /// @{

    // surface grid definition: bin-content: std::array<dindex, 9>
    template <typename grid_shape_t, typename bin_entry_t, typename container_t>
    using surface_grid_t =
        grid<coordinate_axes<grid_shape_t, false, container_t>, bin_entry_t,
             simple_serializer, regular_attacher<9>>;

    // cylindrical grid for the barrel layers
    template <typename bin_entry_t, typename container_t>
    using cylinder_sf_grid =
        surface_grid_t<cylinder2D<>::axes<>, bin_entry_t, container_t>;

    // disc grid for the endcap layers
    template <typename bin_entry_t, typename container_t>
    using disc_sf_grid =
        surface_grid_t<ring2D<>::axes<>, bin_entry_t, container_t>;

    /// @}

    using bfield_backend_t = _bfield_backend_t;

    /// How to store coordinate transform matrices
    template <template <typename...> class vector_t = dvector>
    using transform_store = single_store<__plugin::transform3<detray::scalar>,
                                         vector_t, geometry_context>;

    /// Mask type ids
    enum class mask_ids {
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_portal_cylinder2 = 2,
        e_portal_ring2 = 3,
        e_cylinder2 = 2,
    };

    /// How to store masks
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using mask_store =
        regular_multi_store<mask_ids, empty_context, tuple_t, vector_t,
                            rectangle, trapezoid, cylinder_portal, disc_portal>;

    /// Material type ids
    enum class material_ids {
        e_slab = 0,
        e_none = 1,
    };

    /// How to store materials
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using material_store = regular_multi_store<material_ids, empty_context,
                                               tuple_t, vector_t, slab>;

    /// How to link to the entries in the data stores
    using transform_link = typename transform_store<>::link_type;
    using mask_link = typename mask_store<>::single_link;
    using material_link = typename material_store<>::single_link;
    using source_link = dindex;
    /// Surface type used for sensitives, passives and portals
    using surface_type =
        surface_descriptor<mask_link, material_link, transform_link, nav_link,
                           source_link>;

    /// Portals and passives in the brute froce search, sensitives in the grids
    enum geo_objects : std::size_t {
        e_sensitive = 1,
        e_portal = 0,
        e_passive = 0,
        e_size = 2,
        e_all = e_size,
    };

    /// Acceleration data structures
    enum class sf_finder_ids {
        e_brute_force = 0,     // test all surfaces in a volume (brute force)
        e_disc_grid = 1,       // endcap
        e_cylinder2_grid = 2,  // barrel
        e_default = e_brute_force,
    };

    /// One link for portals/passives and one sensitive surfaces
    using object_link_type =
        dmulti_index<dtyped_index<sf_finder_ids, dindex>, geo_objects::e_size>;

    /// How to store the acceleration data structures
    template <template <typename...> class tuple_t = dtuple,
              typename container_t = host_container_types>
    using surface_finder_store = multi_store<
        sf_finder_ids, empty_context, tuple_t,
        brute_force_collection<surface_type, container_t>,
        grid_collection<disc_sf_grid<surface_type, container_t>>,
        grid_collection<cylinder_sf_grid<surface_type, container_t>>>;

    /// Volume search grid
    template <typename container_t = host_container_types>
    using volume_finder =
        grid<coordinate_axes<
                 cylinder3D::axes<n_axis::bounds::e_open, n_axis::irregular,
                                  n_axis::regular, n_axis::irregular>,
                 true, container_t>,
             dindex, simple_serializer, replacer>;
};

}  // namespace detray
