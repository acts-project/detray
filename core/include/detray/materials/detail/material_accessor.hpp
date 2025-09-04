/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/type_traits.hpp"

namespace detray::detail {

/// @brief Access a single unit of material in different types of material
/// description
namespace material_accessor {

/// Access to material slabs or rods in a homogeneous material
/// description and to raw material in a homogeneous volume material
/// description
template <class material_coll_t, typename point_t = void,
          std::enable_if_t<
              detail::is_hom_material<typename material_coll_t::value_type>::value ,
              bool> = true>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    const material_coll_t &material_coll, const dindex idx,
    const point_t &) noexcept {

    return material_coll[idx];
}

/// Access to material slabs in a material map or volume material
template <
    class material_coll_t,
    std::enable_if_t<detail::is_grid<typename material_coll_t::value_type>::value ,
                     bool> = true>
DETRAY_HOST_DEVICE inline constexpr decltype(auto) get(
    const material_coll_t &material_coll, const dindex idx,
    const typename material_coll_t::value_type::point_type
        &loc_point) noexcept {

    // Find the material slab (only one entry per bin)
    return *(material_coll[idx].search(loc_point));
}

}  // namespace material_accessor

}  // namespace detray::detail
