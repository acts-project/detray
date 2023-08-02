/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/cylinder_portal_intersector.hpp"
#include "detray/io/common/detail/definitions.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/utils/tuple_helpers.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace detray::detail {

/// Determine the type and id of a shape of a mask without triggering a compiler
/// error (sfinae) if the detector does not know the type / enum entry
/// @{
/// Mask shape unknown by detector
template <io::detail::mask_shape shape, typename detector_t, typename = void>
struct mask_info {
    using type = void;
    static constexpr int value{-1};
};

/// Check for a stereo annulus shape
template <typename detector_t>
struct mask_info<io::detail::mask_shape::annulus2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<annulus2D<>, std::uint_least16_t>>(),
                                  void>> {
    using type = annulus2D<>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_annulus2};
};

/// Check for a 2D cylinder shape
template <typename detector_t>
struct mask_info<io::detail::mask_shape::cylinder2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<mask<
                                      cylinder2D<>, std::uint_least16_t>>(),
                                  void>> {
    using type = cylinder2D<>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_cylinder2};
};

/// Check for a 2D cylinder portal shape
template <typename detector_t>
struct mask_info<
    io::detail::mask_shape::portal_cylinder2, detector_t,
    std::enable_if_t<detector_t::masks::template is_defined<
                         mask<cylinder2D<false, cylinder_portal_intersector>,
                              std::uint_least16_t>>(),
                     void>> {
    using type = cylinder2D<false, cylinder_portal_intersector>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_portal_cylinder2};
};

/// Check for a cell wire line shape
template <typename detector_t>
struct mask_info<io::detail::mask_shape::cell_wire, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<line<true>, std::uint_least16_t>>(),
                                  void>> {
    using type = line<true>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_cell_wire};
};

/// Check for a straw wire line shape
template <typename detector_t>
struct mask_info<io::detail::mask_shape::straw_wire, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<line<false>, std::uint_least16_t>>(),
                                  void>> {
    using type = line<false>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_straw_wire};
};

/// Check for a rectangle shape
template <typename detector_t>
struct mask_info<io::detail::mask_shape::rectangle2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<mask<
                                      rectangle2D<>, std::uint_least16_t>>(),
                                  void>> {
    using type = rectangle2D<>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_rectangle2};
};

/// Check for a ring/disc shape
template <typename detector_t>
struct mask_info<io::detail::mask_shape::ring2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<ring2D<>, std::uint_least16_t>>(),
                                  void>> {
    using type = ring2D<>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_portal_ring2};
};

/// Check for a single masked value (1st value is checked)
template <typename detector_t>
struct mask_info<io::detail::mask_shape::single1, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<single3D<0>, std::uint_least16_t>>(),
                                  void>> {
    using type = single3D<0>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_single1};
};

/// Check for a single masked value (2nd value is checked)
template <typename detector_t>
struct mask_info<io::detail::mask_shape::single2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<single3D<1>, std::uint_least16_t>>(),
                                  void>> {
    using type = single3D<1>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_single2};
};

/// Check for a single masked value (3rd value is checked)
template <typename detector_t>
struct mask_info<io::detail::mask_shape::single3, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<single3D<2>, std::uint_least16_t>>(),
                                  void>> {
    using type = single3D<2>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_single3};
};

/// Check for a trapezoid shape
template <typename detector_t>
struct mask_info<io::detail::mask_shape::trapezoid2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<mask<
                                      trapezoid2D<>, std::uint_least16_t>>(),
                                  void>> {
    using type = trapezoid2D<>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_trapezoid2};
};
/// @}

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

/// Is the value type in the detector material store a simple material or is it
/// wrapped in another class (e.g. grids for material maps)
template <class detector_t>
struct has_homogeneous_material<
    detector_t, std::enable_if_t<has_material_slabs_v<detector_t> or
                                     has_material_rods_v<detector_t>,
                                 void>> : public std::true_type {};

template <typename T>
inline constexpr bool has_homogeneous_material_v =
    has_homogeneous_material<T>::value;

}  // namespace detray::detail
