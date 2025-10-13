/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/core/concepts.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/tracking_volume.hpp"

// System include(s).
#include <vector>

namespace detray {

/// @returns the total number of portals in the detector @param det
template <typename detector_t>
DETRAY_HOST_DEVICE inline std::size_t n_portals(const detector_t& det) {
    std::size_t n_portals{0u};

    for (const auto& vol_desc : det.volumes()) {
        n_portals += tracking_volume{det, vol_desc}.n_portals();
    }

    return n_portals;
}

/// @returns the total number of sensitive surfaces in the detector @param det
template <typename detector_t>
DETRAY_HOST_DEVICE inline std::size_t n_sensitives(const detector_t& det) {
    std::size_t n_sensitives{0u};

    for (const auto& vol_desc : det.volumes()) {
        n_sensitives += tracking_volume{det, vol_desc}.n_sensitives();
    }

    return n_sensitives;
}

/// @returns the total number of passive surfaces in the detector @param det
template <typename detector_t>
DETRAY_HOST_DEVICE inline std::size_t n_passives(const detector_t& det) {
    std::size_t n_passives{0u};

    for (const auto& vol_desc : det.volumes()) {
        n_passives += tracking_volume{det, vol_desc}.n_passives();
    }

    return n_passives;
}

/// @returns the total number of grids in a detector store @param store
template <std::size_t I = 0u, typename store_t>
DETRAY_HOST_DEVICE inline std::size_t n_grids(const store_t& store,
                                              std::size_t n = 0u) {

    constexpr auto coll_id{store_t::value_types::to_id(I)};
    using value_type = typename store_t::template get_type<coll_id>;

    if constexpr (detray::concepts::grid<value_type>) {
        n += store.template size<coll_id>();
    }

    if constexpr (I < store_t::n_collections() - 1u) {
        return n_grids<I + 1>(store, n);
    }
    return n;
}

/// @returns the total number of surface grids in the detector  @param det
template <typename detector_t>
DETRAY_HOST_DEVICE inline std::size_t n_surface_grids(const detector_t& det) {
    return n_grids(det.accelerator_store());
}

/// @returns the total number of material maps in the detector  @param det
template <typename detector_t>
DETRAY_HOST_DEVICE inline std::size_t n_material_maps(const detector_t& det) {
    return n_grids(det.material_store());
}

/// @returns the total number of material slabs outside of material maps in the
/// detector  @param det
template <typename detector_t>
DETRAY_HOST_DEVICE inline std::size_t n_material_slabs(const detector_t& det) {
    if constexpr (detray::concepts::has_material_slabs<detector_t>) {
        constexpr auto slab_id{detector_t::materials::id::e_slab};
        return det.material_store().template size<slab_id>();
    } else {
        return 0u;
    }
}

/// @returns the total number of material rods in the detector  @param det
template <typename detector_t>
DETRAY_HOST_DEVICE inline std::size_t n_material_rods(const detector_t& det) {
    if constexpr (detray::concepts::has_material_rods<detector_t>) {
        constexpr auto rod_id{detector_t::materials::id::e_rod};
        return det.material_store().template size<rod_id>();
    } else {
        return 0u;
    }
}

}  // namespace detray
