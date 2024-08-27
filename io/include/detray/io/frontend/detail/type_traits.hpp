/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace detray::detail {

/// Check for grid types in a data store
/// @{
template <typename>
struct contains_grids;

template <template <typename, typename...> class R, class ID, typename... Ts>
struct contains_grids<R<ID, Ts...>> {
    static constexpr bool value{std::disjunction_v<detail::is_grid<Ts>...>};
};

template <typename T>
inline constexpr bool contains_grids_v = contains_grids<T>::value;
/// @}

/// Check for the various types of material
/// @{
template <class detector_t, typename = void>
struct has_surface_grids : public std::false_type {};

template <class detector_t>
struct has_surface_grids<
    detector_t,
    std::enable_if_t<contains_grids_v<typename detector_t::accel>, void>>
    : public std::true_type {};

template <typename T>
inline constexpr bool has_surface_grids_v = has_surface_grids<T>::value;
/// @}

/// Check for the various types of material
/// @{
template <class detector_t, typename = void>
struct has_material_slabs : public std::false_type {};

template <class detector_t>
struct has_material_slabs<
    detector_t,
    std::enable_if_t<detector_t::materials::template is_defined<
                         material_slab<typename detector_t::scalar_type>>(),
                     void>> : public std::true_type {};

template <typename T>
inline constexpr bool has_material_slabs_v = has_material_slabs<T>::value;

template <class detector_t, typename = void>
struct has_material_rods : public std::false_type {};

template <class detector_t>
struct has_material_rods<
    detector_t,
    std::enable_if_t<detector_t::materials::template is_defined<
                         material_rod<typename detector_t::scalar_type>>(),
                     void>> : public std::true_type {};

template <typename T>
inline constexpr bool has_material_rods_v = has_material_rods<T>::value;

template <class detector_t, typename = void>
struct has_homogeneous_material : public std::false_type {};

template <class detector_t>
struct has_homogeneous_material<
    detector_t, std::enable_if_t<has_material_slabs_v<detector_t> ||
                                     has_material_rods_v<detector_t>,
                                 void>> : public std::true_type {};

template <typename T>
inline constexpr bool has_homogeneous_material_v =
    has_homogeneous_material<T>::value;

template <class detector_t, typename = void>
struct has_material_grids : public std::false_type {};

template <class detector_t>
struct has_material_grids<
    detector_t,
    std::enable_if_t<contains_grids_v<typename detector_t::materials>, void>>
    : public std::true_type {};

template <typename T>
inline constexpr bool has_material_grids_v = has_material_grids<T>::value;
/// @}

}  // namespace detray::detail
