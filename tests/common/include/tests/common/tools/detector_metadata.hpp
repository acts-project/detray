/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/core/surfaces_finder.hpp"
#include "detray/core/transform_store.hpp"
#include "detray/core/type_registry.hpp"
#include "detray/grids/axis.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"
#include "detray/masks/masks.hpp"
#include "detray/utils/indexing.hpp"

namespace detray {

struct volume_stats {
    std::size_t n_max_objects_per_volume = 0;
};

/// edge links: next volume, next (local) object finder
using edge_type = std::array<dindex, 2>;

/// mask types
using rectangle = rectangle2<planar_intersector,
                             __plugin::cartesian2<detray::scalar>, edge_type>;
using trapezoid = trapezoid2<planar_intersector,
                             __plugin::cartesian2<detray::scalar>, edge_type>;
using annulus = annulus2<planar_intersector,
                         __plugin::cartesian2<detray::scalar>, edge_type>;
using cylinder = cylinder3<false, cylinder_intersector,
                           __plugin::cylindrical2<detray::scalar>, edge_type>;
using disc =
    ring2<planar_intersector, __plugin::cartesian2<detray::scalar>, edge_type>;

/// Defines all available types
template <typename dynamic_data, std::size_t NGRIDS = 1>
struct full_metadata {

    // How many grids have to be built
    enum grids : std::size_t {
        n_grids = NGRIDS,
    };

    // How to store and link transforms
    template <template <typename...> class vector_t = dvector>
    using transform_store = static_transform_store<vector_t>;

    template <typename... objects>
    using object_definitions = object_registry<objects...>;

    /// Give your mask types a name (needs to be consecutive to be matched
    /// to a type!)
    enum mask_ids : unsigned int {
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_annulus2 = 2,
        e_cylinder3 = 3,
        e_portal_cylinder3 = 3,  // no distinction from surface cylinder
        e_ring2 = 4,
        e_portal_ring2 = 4,
        e_single3 = 5,
    };

    // How to store and link masks
    using mask_definitions =
        mask_registry<mask_ids, rectangle, trapezoid, annulus, cylinder, disc>;

    // Accelerator types
    template <template <typename, std::size_t> class array_t = darray,
              template <typename...> class vector_t = dvector,
              template <typename...> class tuple_t = dtuple,
              template <typename...> class jagged_vector_t = djagged_vector>
    using volume_finder =
        grid2<replace_populator, axis::irregular, axis::irregular, serializer2,
              vector_t, jagged_vector_t, array_t, tuple_t, dindex>;

    template <template <typename, std::size_t> class array_t = darray,
              template <typename...> class vector_t = dvector,
              template <typename...> class tuple_t = dtuple,
              template <typename...> class jagged_vector_t = djagged_vector>
    using surface_finder =
        surfaces_finder<n_grids, array_t, tuple_t, vector_t, jagged_vector_t>;

    dynamic_data _data;
};

/// Defines the data types needed for the toy detector
struct toy_metadata {

    // How many grids have to be built
    enum grids : std::size_t {
        n_grids = 20,
    };

    // How to store and link transforms
    template <template <typename...> class vector_t = dvector>
    using transform_store = static_transform_store<vector_t>;

    template <typename... objects>
    using object_definitions = object_registry<objects...>;

    /// Give your mask types a name (needs to be consecutive to be matched
    /// to a type!)
    enum mask_ids : unsigned int {
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_cylinder3 = 2,         // Put the beampipe into the same container as
        e_portal_cylinder3 = 2,  // the cylinder portals
        e_portal_ring2 = 3,
    };

    // How to store and link masks
    using mask_definitions =
        mask_registry<mask_ids, rectangle, trapezoid, cylinder, disc>;

    // Accelerator types
    template <template <typename, std::size_t> class array_t = darray,
              template <typename...> class vector_t = dvector,
              template <typename...> class tuple_t = dtuple,
              template <typename...> class jagged_vector_t = djagged_vector>
    using volume_finder =
        grid2<replace_populator, axis::irregular, axis::irregular, serializer2,
              vector_t, jagged_vector_t, array_t, tuple_t, dindex>;

    template <template <typename, std::size_t> class array_t = darray,
              template <typename...> class vector_t = dvector,
              template <typename...> class tuple_t = dtuple,
              template <typename...> class jagged_vector_t = djagged_vector>
    using surface_finder =
        surfaces_finder<n_grids, array_t, tuple_t, vector_t, jagged_vector_t>;

    volume_stats _data;
};

/// Defines a detector with only rectangle surfaces
struct telescope_metadata {

    // How many grids have to be built
    enum grids : std::size_t {
        n_grids = 1,
    };

    // How to store and link transforms
    template <template <typename...> class vector_t = dvector>
    using transform_store = static_transform_store<vector_t>;

    template <typename... objects>
    using object_definitions = object_registry<objects...>;

    /// Give your mask types a name (needs to be consecutive to be matched
    /// to a type!)
    enum mask_ids : unsigned int {
        e_rectangle2 = 0,
    };

    // How to store and link masks
    using mask_definitions = mask_registry<mask_ids, rectangle>;

    // Accelerator types
    using volume_finder = void;
    using surface_finder = void;
};

struct detector_registry {
    using default_detector = full_metadata<volume_stats, 1>;
    using tml_detector = full_metadata<volume_stats, 192>;
    using toy_detector = toy_metadata;
    using telescope_detector = telescope_metadata;
};

}  // namespace detray