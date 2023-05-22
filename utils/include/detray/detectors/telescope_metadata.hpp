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
#include "detray/intersection/plane_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/surface_finders/brute_force_finder.hpp"

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/vector.hpp>

namespace detray {

/// Defines a detector with only rectangle/unbounded surfaces
template <typename mask_shape_t = rectangle2D<>,
          typename _bfield_backend_t =
              covfie::backend::constant<covfie::vector::vector_d<scalar, 3>,
                                        covfie::vector::vector_d<scalar, 3>>>
struct telescope_metadata {

    /// mask to (next) volume link: next volume(s)
    using nav_link = std::uint_least16_t;

    /// mask types (these types are needed for the portals, which are always
    /// there, and to resolve the material, i.e. slab vs. rod)
    using rectangle = mask<rectangle2D<>, nav_link>;
    using straw_wire = mask<line<false>, nav_link>;
    using cell_wire = mask<line<true>, nav_link>;

    /// material types
    using rod = material_rod<detray::scalar>;
    using slab = material_slab<detray::scalar>;

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
        e_rectangle2 = 0,
        e_portal_rectangle2 = 0,
        e_annulus2 = 1,
        e_cell_wire = 1,
        e_cylinder2 = 1,
        e_ring2 = 1,
        e_trapezoid2 = 1,
        e_single1 = 1,
        e_single2 = 1,
        e_single3 = 1,
        e_straw_wire = 1,
        e_unbounded_annulus2 = 1,
        e_unbounded_cell2 = 1,
        e_unbounded_cylinder2 = 1,
        e_unbounded_disc2 = 1,
        e_unbounded_rectangle2 = 1,
        e_unbounded_trapezoid2 = 1,
        e_unbounded_straw2 = 1,
        e_unmasked2 = 1,
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
        std::conditional_t<std::is_same_v<mask<mask_shape_t>, cell_wire> |
                               std::is_same_v<mask<mask_shape_t>, straw_wire>,
                           regular_multi_store<material_ids, empty_context,
                                               tuple_t, vector_t, slab, rod>,
                           regular_multi_store<material_ids, empty_context,
                                               tuple_t, vector_t, slab>>;

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
    using volume_finder = brute_force_collection<surface_type, container_t>;
};

}  // namespace detray
