/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
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
#include "detray/intersection/plane_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/unmasked.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/surface_finders/accelerator_grid.hpp"
#include "detray/surface_finders/brute_force_finder.hpp"

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/vector.hpp>

namespace detray {

struct volume_stats {
    std::size_t n_max_objects_per_volume = 0;
};

/// mask to (next) volume link: next volume(s)
using volume_link_type = dindex;

/// mask types
using rectangle = mask<rectangle2D<>, volume_link_type>;
using trapezoid = mask<trapezoid2D<>, volume_link_type>;
using annulus = mask<annulus2D<>, volume_link_type>;
using cylinder = mask<cylinder2D<>, volume_link_type>;
using disc = mask<ring2D<>, volume_link_type>;
using unbounded_plane = mask<unmasked, volume_link_type>;
using lines = mask<line<>, volume_link_type>;

/// material types
using slab = material_slab<detray::scalar>;
using rod = material_rod<detray::scalar>;

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
using disc_sf_grid = surface_grid_t<ring2D<>::axes<>, bin_entry_t, container_t>;

/// @}

/// Defines all available types
template <typename dynamic_data, std::size_t kBrlGrids = 1,
          std::size_t kEdcGrids = 1, std::size_t kDefault = 1,
          typename _bfield_backend_t =
              covfie::backend::constant<covfie::vector::vector_d<scalar, 3>,
                                        covfie::vector::vector_d<scalar, 3>>>
struct full_metadata {
    using bfield_backend_t = _bfield_backend_t;

    /// How to index the constituent objects in a volume
    /// If they share the same index value here, they will be added into the
    /// same container range without any sorting guarantees
    enum geo_objects : std::size_t {
        e_sensitive = 0,
        e_portal = 0,
        e_passive = 0,
        e_size = 1,
        e_all = e_size,
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
        e_trapezoid2 = 1,
        e_annulus2 = 2,
        e_cylinder2 = 3,
        e_portal_cylinder2 = 3,  // no distinction from surface cylinder
        e_ring2 = 4,
        e_portal_ring2 = 4,
    };

    /// How to store and link masks
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using mask_store =
        regular_multi_store<mask_ids, empty_context, tuple_t, vector_t,
                            rectangle, trapezoid, annulus, cylinder, disc>;

    /// Give your material types a name (needs to be consecutive to be matched
    /// to a type!)
    enum class material_ids {
        e_slab = 0,
        e_rod = 1,
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
    using surface_type =
        surface<mask_link, material_link, transform_link, source_link>;

    /// Surface finders
    enum class sf_finder_ids {
        e_brute_force = 0,  // test all surfaces in a volume (brute force)
        // e_disc_grid = 1,      // barrel
        // e_cylinder_grid = 2,  // endcap
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
  grid_collection<cylinder_sf_grid<surface_type,
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

/// Defines the data types needed for the toy detector
template <typename _bfield_backend_t =
              covfie::backend::constant<covfie::vector::vector_d<scalar, 3>,
                                        covfie::vector::vector_d<scalar, 3>>>
struct toy_metadata {
    using bfield_backend_t = _bfield_backend_t;

    /// How to index the constituent objects in a volume
    /// If they share the same index value here, they will be added into the
    /// same container range without any sorting guarantees
    enum geo_objects : std::size_t {
        e_sensitive = 0,
        e_portal = 0,
        e_passive = 0,
        e_size = 1,
        e_all = e_size,
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
        e_trapezoid2 = 1,
        e_cylinder2 = 2,         // Put the beampipe into the same container as
        e_portal_cylinder2 = 2,  // the cylinder portals
        e_portal_ring2 = 3,
    };

    /// How to store and link masks
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using mask_store =
        regular_multi_store<mask_ids, empty_context, tuple_t, vector_t,
                            rectangle, trapezoid, cylinder, disc>;

    /// Give your material types a name (needs to be consecutive to be matched
    /// to a type!)
    enum class material_ids {
        e_slab = 0,
        e_none = 1,
    };

    /// How to store and link materials
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using material_store = regular_multi_store<material_ids, empty_context,
                                               tuple_t, vector_t, slab>;

    /// Surface type used for sensitives, passives and portals
    using transform_link = typename transform_store<>::link_type;
    using mask_link = typename mask_store<>::single_link;
    using material_link = typename material_store<>::single_link;
    using source_link = dindex;
    using surface_type =
        surface<mask_link, material_link, transform_link, source_link>;

    /// Surface finders
    enum class sf_finder_ids {
        e_brute_force = 0,    // test all surfaces in a volume (brute force)
        e_disc_grid = 1,      // barrel
        e_cylinder_grid = 2,  // endcap
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
  grid_collection<cylinder_sf_grid<surface_type,
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

/// Defines a detector with only rectangle/unbounded surfaces
template <typename mask_shape_t = rectangle2D<>,
          typename _bfield_backend_t =
              covfie::backend::constant<covfie::vector::vector_d<scalar, 3>,
                                        covfie::vector::vector_d<scalar, 3>>>
struct telescope_metadata {
    using bfield_backend_t = _bfield_backend_t;

    /// How to index the constituent objects in a volume
    /// If they share the same index value here, they will be added into the
    /// same container range without any sorting guarantees
    enum geo_objects : std::size_t {
        e_sensitive = 0,
        e_portal = 0,
        e_size = 1,
        e_all = e_size,
    };

    /// How a volume finds its constituent objects in the detector containers
    /// In this case: One range for sensitive/passive surfaces, one for portals
    using object_link_type = dmulti_index<dindex_range, geo_objects::e_size>;

    /// How to store and link transforms
    template <template <typename...> class vector_t = dvector>
    using transform_store = single_store<__plugin::transform3<detray::scalar>,
                                         vector_t, geometry_context>;

    /// Rectangles are always needed as portals. Only one mask shape is allowed
    /// in addition
    enum class mask_ids {
        e_portal_rectangle2 = 0,
        e_rectangle2 = 0,
        e_annulus2 = 1,
        e_cylinder2 = 1,
        e_ring2 = 1,
        e_trapezoid2 = 1,
        e_unbounded = 1
    };

    /// How to store and link masks
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using mask_store = std::conditional_t<
        std::is_same_v<mask<mask_shape_t>, rectangle>,
        regular_multi_store<mask_ids, empty_context, tuple_t, vector_t,
                            rectangle>,
        regular_multi_store<mask_ids, empty_context, tuple_t, vector_t,
                            rectangle, mask<mask_shape_t>>>;

    /// Give your material types a name (needs to be consecutive to be matched
    /// to a type!)
    enum class material_ids {
        e_slab = 0,
        e_rod = 0,
        e_none = 1,
    };

    /// How to store and link materials
    template <template <typename...> class tuple_t = dtuple,
              template <typename...> class vector_t = dvector>
    using material_store =
        std::conditional_t<std::is_same_v<mask<mask_shape_t>, lines>,
                           regular_multi_store<material_ids, empty_context,
                                               tuple_t, vector_t, slab, rod>,
                           regular_multi_store<material_ids, empty_context,
                                               tuple_t, vector_t, slab>>;

    /// Surface type used for sensitives, passives and portals
    using transform_link = typename transform_store<>::link_type;
    using mask_link = typename mask_store<>::single_link;
    using material_link = typename material_store<>::single_link;
    using source_link = dindex;
    using surface_type =
        surface<mask_link, material_link, transform_link, source_link>;

    /// Surface finders
    enum class sf_finder_ids {
        e_brute_force = 0,  // test all surfaces in a volume (brute force)
        e_default = e_brute_force,
    };

    /// How to store and link surface grids
    template <template <typename...> class tuple_t = dtuple,
              typename container_t = host_container_types>
    using surface_finder_store =
        multi_store<sf_finder_ids, empty_context, tuple_t,
                    brute_force_collection<surface_type, container_t>>;

    /// Volume grid
    template <typename container_t = host_container_types>
    using volume_finder =
        grid<coordinate_axes<
                 cylinder3D::axes<n_axis::bounds::e_open, n_axis::irregular,
                                  n_axis::regular, n_axis::irregular>,
                 true, container_t>,
             dindex, simple_serializer, replacer>;
};

struct detector_registry {
    using default_detector = full_metadata<volume_stats, 1>;
    using tml_detector = full_metadata<volume_stats, 192>;
    using toy_detector = toy_metadata<>;
    template <typename mask_shape_t>
    using telescope_detector = telescope_metadata<mask_shape_t>;
};

}  // namespace detray