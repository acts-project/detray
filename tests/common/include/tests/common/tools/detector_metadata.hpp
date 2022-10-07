/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/surfaces_finder.hpp"
#include "detray/core/transform_store.hpp"
#include "detray/core/type_registry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/grids/axis.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"

namespace detray {

struct volume_stats {
    std::size_t n_max_objects_per_volume = 0;
};

/// edge links: next volume, next (local) object finder
using edge_type = std::array<dindex, 2>;

/// mask types
using rectangle = rectangle2<__plugin::transform3<detray::scalar>,
                             plane_intersector, cartesian2, edge_type>;
using trapezoid = trapezoid2<__plugin::transform3<detray::scalar>,
                             plane_intersector, cartesian2, edge_type>;
// @TODO: Use Polar2?
using annulus = annulus2<__plugin::transform3<detray::scalar>,
                         plane_intersector, polar2, edge_type>;
using cylinder = cylinder3<__plugin::transform3<detray::scalar>,
                           cylinder_intersector, cylindrical2, edge_type>;
// @TODO: Use Polar2?
using disc = ring2<__plugin::transform3<detray::scalar>, plane_intersector,
                   polar2, edge_type>;
using unbounded_plane = unmasked<__plugin::transform3<detray::scalar>,
                                 plane_intersector, cartesian2, edge_type>;

using slab = material_slab<detray::scalar>;
using rod = material_rod<detray::scalar>;

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
        tuple_vector_registry<mask_ids, rectangle, trapezoid, annulus, cylinder,
                              disc>;

    /// Give your material types a name (needs to be consecutive to be matched
    /// to a type!)
    enum material_ids : unsigned int {
        e_slab = 0,
        e_rod = 1,
    };

    // How to store and link materials
    using material_definitions = tuple_vector_registry<material_ids, slab, rod>;

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
        tuple_vector_registry<mask_ids, rectangle, trapezoid, cylinder, disc>;

    /// Give your material types a name (needs to be consecutive to be matched
    /// to a type!)
    enum material_ids : unsigned int {
        e_slab = 0,
    };

    // How to store and link materials
    using material_definitions = tuple_vector_registry<material_ids, slab>;

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

/// Defines a detector with only rectangle/unbounded surfaces
struct telescope_metadata {

    // How many grids have to be built
    enum grids : std::size_t {
        n_grids = 0,
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
        e_unbounded_plane2 = 1,
    };

    // How to store and link masks
    using mask_definitions =
        tuple_vector_registry<mask_ids, rectangle, unbounded_plane>;

    /// Give your material types a name (needs to be consecutive to be matched
    /// to a type!)
    enum material_ids : unsigned int {
        e_slab = 0,
    };

    // How to store and link materials
    using material_definitions = tuple_vector_registry<material_ids, slab>;

    // Accelerator types (are not used)
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
};

struct detector_registry {
    using default_detector = full_metadata<volume_stats, 1>;
    using tml_detector = full_metadata<volume_stats, 192>;
    using toy_detector = toy_metadata;
    using telescope_detector = telescope_metadata;
};

}  // namespace detray