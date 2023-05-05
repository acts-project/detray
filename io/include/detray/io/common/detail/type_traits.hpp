/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/utils/tuple_helpers.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace detray::detail {

template <class detector_t, typename = void>
struct is_homogeneous_material : public std::false_type {};

/// Is the value type in the detector material store a simple material or is it
/// wrapped in another class (e.g. grids for material maps)
template <class detector_t>
struct is_homogeneous_material<
    detector_t,
    std::enable_if_t<
        std::is_base_of_v<detail::homogeneous_material_tag,
                          typename detail::tuple_element_t<
                              0, typename detector_t::material_container::
                                     tuple_type>::value_type>,
        void>> : public std::true_type {};

template <typename T>
inline constexpr bool is_homogeneous_material_v =
    is_homogeneous_material<T>::value;

}  // namespace detray::detail
