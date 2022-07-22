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
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/surface_finders/brute_force_finder.hpp"
#include "detray/surface_finders/grid2_finder.hpp"

namespace detray {

struct volume_stats {
    std::size_t n_max_objects_per_volume = 0;
};

/// edge links: next volume, next (local) object finder
using edge_type = std::array<dindex, 2>;

/// mask types
using rectangle = rectangle2<__plugin::cartesian2<detray::scalar>, edge_type>;
using trapezoid = trapezoid2<__plugin::cartesian2<detray::scalar>, edge_type>;
using annulus = annulus2<__plugin::cartesian2<detray::scalar>, edge_type>;
using cylinder = cylinder3<cylinder_intersector,
                           __plugin::cylindrical2<detray::scalar>, edge_type>;
using disc = ring2<__plugin::cartesian2<detray::scalar>, edge_type>;
using unbounded_plane =
    unmasked<__plugin::cartesian2<detray::scalar>, edge_type>;

using slab = material_slab<detray::scalar>;
using rod = material_rod<detray::scalar>;

/// Defines all available types
template <typename dynamic_data, std::size_t kBrlGrids = 1,
          std::size_t kEdcGrids = 1, std::size_t kDefault = 1>
struct full_metadata {

    /// How to store and link transforms
    template <template <typename...> class vector_t = dvector>
    using transform_store = static_transform_store<vector_t>;

    template <typename... objects>
    using object_definitions = object_registry<objects...>;

    /// Give your mask types a name (needs to be consecutive to be matched
    /// to a type!)
    enum class mask_ids {
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_annulus2 = 2,
        e_cylinder3 = 3,
        e_portal_cylinder3 = 3,  // no distinction from surface cylinder
        e_ring2 = 4,
        e_portal_ring2 = 4,
        e_single3 = 5,
    };

    /// How to store and link masks
    using mask_definitions =
        tuple_vector_registry<mask_ids, rectangle, trapezoid, annulus, cylinder,
                              disc>;

    /// Give your material types a name (needs to be consecutive to be matched
    /// to a type!)
    enum class material_ids {
        e_slab = 0,
        e_rod = 1,
    };

    // How to store and link materials
    using material_definitions = tuple_vector_registry<material_ids, slab, rod>;

    /// How many grids have to be built
    enum grids : std::size_t {
        n_other = kDefault,
        n_z_phi_grids = kBrlGrids,
        n_r_phi_grids = kEdcGrids
    };

    /// Surface finders
    enum class sf_finder_ids {
        e_brute_force = 0,  // test all surfaces in a volume (brute force)
        e_z_phi_grid = 1,   // barrel
        e_r_phi_grid = 2,   // endcap
        e_default = e_brute_force,
    };

    /// Surface finder types
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
        surfaces_finder<n_z_phi_grids + n_r_phi_grids, array_t, tuple_t,
                        vector_t, jagged_vector_t>;

    template <template <typename, std::size_t> class array_t = darray,
              template <typename...> class vector_t = dvector,
              template <typename...> class tuple_t = dtuple,
              template <typename...> class jagged_vector_t = djagged_vector>
    using sf_finder_definitions = tuple_array_registry<
        sf_finder_ids,
        std::index_sequence<n_other, n_z_phi_grids, n_r_phi_grids>,
        brute_force_finder,
        regular_circular_grid<vector_t, jagged_vector_t, array_t, tuple_t>,
        regular_circular_grid2<vector_t, jagged_vector_t, array_t, tuple_t>>;

    dynamic_data _data;
};

/// Defines the data types needed for the toy detector
struct toy_metadata {

    /// How to store and link transforms
    template <template <typename...> class vector_t = dvector>
    using transform_store = static_transform_store<vector_t>;

    template <typename... objects>
    using object_definitions = object_registry<objects...>;

    /// Give your mask types a name (needs to be consecutive to be matched
    /// to a type!)
    enum class mask_ids {
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_cylinder3 = 2,         // Put the beampipe into the same container as
        e_portal_cylinder3 = 2,  // the cylinder portals
        e_portal_ring2 = 3,
    };

    /// How to store and link masks
    using mask_definitions =
        tuple_vector_registry<mask_ids, rectangle, trapezoid, cylinder, disc>;

    /// Give your material types a name (needs to be consecutive to be matched
    /// to a type!)
    enum class material_ids {
        e_slab = 0,
    };

    /// How to store and link materials
    using material_definitions = tuple_vector_registry<material_ids, slab>;

    /// How many grids have to be built
    enum grids : std::size_t {
        n_other = 1,
        n_barrel_grids = 4,
        n_endcap_grids = 14,
    };

    /// Surface finders
    enum class sf_finder_ids {
        e_brute_force = 0,  // test all surfaces in a volume (brute force)
        e_z_phi_grid = 1,   // barrel
        e_r_phi_grid = 2,   // endcap
        e_default = e_brute_force,
    };

    ///  Surface finder types
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
        surfaces_finder<n_barrel_grids + n_endcap_grids, array_t, tuple_t,
                        vector_t, jagged_vector_t>;

    template <template <typename, std::size_t> class array_t = darray,
              template <typename...> class vector_t = dvector,
              template <typename...> class tuple_t = dtuple,
              template <typename...> class jagged_vector_t = djagged_vector>
    using sf_finder_definitions = tuple_array_registry<
        sf_finder_ids,
        std::index_sequence<n_other, n_barrel_grids, n_endcap_grids>,
        brute_force_finder,
        regular_circular_grid<vector_t, jagged_vector_t, array_t, tuple_t>,
        regular_circular_grid2<vector_t, jagged_vector_t, array_t, tuple_t>>;

    volume_stats _data;
};

/// Defines a detector with only rectangle/unbounded surfaces
struct telescope_metadata {

    /// How to store and link transforms
    template <template <typename...> class vector_t = dvector>
    using transform_store = static_transform_store<vector_t>;

    template <typename... objects>
    using object_definitions = object_registry<objects...>;

    /// Give your mask types a name (needs to be consecutive to be matched
    /// to a type!)
    enum class mask_ids {
        e_rectangle2 = 0,
        e_unbounded_plane2 = 1,
    };

    /// How to store and link masks
    using mask_definitions =
        tuple_vector_registry<mask_ids, rectangle, unbounded_plane>;

    /// Give your material types a name (needs to be consecutive to be matched
    /// to a type!)
    enum class material_ids {
        e_slab = 0,
    };

    /// How to store and link materials
    using material_definitions = tuple_vector_registry<material_ids, slab>;

    /// How many grids have to be built
    enum grids : std::size_t {
        n_other = 1,
    };

    /// Surface finders
    enum class sf_finder_ids {
        e_brute_force = 0,  // test all surfaces in a volume (brute force)
        e_default = e_brute_force,
    };

    ///  Surface finder types (are not used in telescope detector)
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
        surfaces_finder<n_other, array_t, tuple_t, vector_t, jagged_vector_t>;

    template <template <typename, std::size_t> class array_t = darray,
              template <typename...> class vector_t = dvector,
              template <typename...> class tuple_t = dtuple,
              template <typename...> class jagged_vector_t = djagged_vector>
    using sf_finder_definitions =
        tuple_array_registry<sf_finder_ids, std::index_sequence<n_other>,
                             brute_force_finder>;
};

struct detector_registry {
    using default_detector = full_metadata<volume_stats, 1>;
    using tml_detector = full_metadata<volume_stats, 192>;
    using toy_detector = toy_metadata;
    using telescope_detector = telescope_metadata;
};

}  // namespace detray