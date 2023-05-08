/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/cylinder_portal_intersector.hpp"
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
                                      cylinder2D<true>, std::uint_least16_t>>(),
                                  void>> {
    using type = cylinder2D<true>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_cylinder2};
};

/// Check for a 2D cylinder portal shape
template <typename detector_t>
struct mask_info<
    io::detail::mask_shape::cylinder2, detector_t,
    std::enable_if_t<detector_t::masks::template is_defined<
                         mask<cylinder2D<false, cylinder_portal_intersector>,
                              std::uint_least16_t>>(),
                     void>> {
    using type = cylinder2D<true>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_portal_cylinder2};
};

/// Check for a line shape
template <typename detector_t>
struct mask_info<io::detail::mask_shape::line, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<line<>, std::uint_least16_t>>(),
                                  void>> {
    using type = line<>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_line};
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

/// Check for a single masked value
template <typename detector_t>
struct mask_info<io::detail::mask_shape::single3, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<single3D<>, std::uint_least16_t>>(),
                                  void>> {
    using type = single3D<>;
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
