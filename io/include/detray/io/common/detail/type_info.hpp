/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/cylinder_portal_intersector.hpp"
#include "detray/io/frontend/definitions.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/material_map.hpp"
#include "detray/utils/tuple_helpers.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace detray::io::detail {

/// Determine the type and id of a shape of a mask without triggering a compiler
/// error (sfinae) if the detector does not know the type / enum entry
/// @{
/// Mask shape unknown by detector
template <io::shape_id shape, typename detector_t, typename = void>
struct mask_info {
    using shape_id = typename detector_t::masks::id;
    using type = void;
    static constexpr shape_id value{detray::detail::invalid_value<shape_id>()};
};

/// Check for a stereo annulus shape
template <typename detector_t>
struct mask_info<io::shape_id::annulus2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<annulus2D<>, std::uint_least16_t>>(),
                                  void>> {
    using type = annulus2D<>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_annulus2};
};

/// Check for a 2D cylinder shape
template <typename detector_t>
struct mask_info<io::shape_id::cylinder2, detector_t,
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
    io::shape_id::portal_cylinder2, detector_t,
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
struct mask_info<io::shape_id::cell_wire, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<line<true>, std::uint_least16_t>>(),
                                  void>> {
    using type = line<true>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_cell_wire};
};

/// Check for a straw wire line shape
template <typename detector_t>
struct mask_info<io::shape_id::straw_wire, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<line<false>, std::uint_least16_t>>(),
                                  void>> {
    using type = line<false>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_straw_wire};
};

/// Check for a rectangle shape
template <typename detector_t>
struct mask_info<io::shape_id::rectangle2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<mask<
                                      rectangle2D<>, std::uint_least16_t>>(),
                                  void>> {
    using type = rectangle2D<>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_rectangle2};
};

/// Check for a ring/disc shape
template <typename detector_t>
struct mask_info<io::shape_id::ring2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<ring2D<>, std::uint_least16_t>>(),
                                  void>> {
    using type = ring2D<>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_portal_ring2};
};

/// Check for a single masked value (1st value is checked)
template <typename detector_t>
struct mask_info<io::shape_id::single1, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<single3D<0>, std::uint_least16_t>>(),
                                  void>> {
    using type = single3D<0>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_single1};
};

/// Check for a single masked value (2nd value is checked)
template <typename detector_t>
struct mask_info<io::shape_id::single2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<single3D<1>, std::uint_least16_t>>(),
                                  void>> {
    using type = single3D<1>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_single2};
};

/// Check for a single masked value (3rd value is checked)
template <typename detector_t>
struct mask_info<io::shape_id::single3, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<single3D<2>, std::uint_least16_t>>(),
                                  void>> {
    using type = single3D<2>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_single3};
};

/// Check for a trapezoid shape
template <typename detector_t>
struct mask_info<io::shape_id::trapezoid2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<mask<
                                      trapezoid2D<>, std::uint_least16_t>>(),
                                  void>> {
    using type = trapezoid2D<>;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_trapezoid2};
};
/// @}

/// Determine the type and id of a material map without triggering a compiler
/// error (sfinae) if the detector does not know the type / enum entry
/// @{
/// Material map unknown by detector
template <io::material_id mat, typename detector_t, typename = void>
struct mat_map_info {
    using material_id = typename detector_t::materials::id;
    using type = void;
    static constexpr material_id value{
        detray::detail::invalid_value<material_id>()};
};

/// Check for a 2D disc material map
template <typename detector_t>
struct mat_map_info<
    io::material_id::ring2_map, detector_t,
    std::enable_if_t<detector_t::materials::template is_defined<material_map<
                         ring2D<>, typename detector_t::scalar_type>>(),
                     void>> {
    using type = material_map<ring2D<>, typename detector_t::scalar_type>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_disc2_map};
};

/// Check for a 2D cartesian material map
template <typename detector_t>
struct mat_map_info<
    io::material_id::rectangle2_map, detector_t,
    std::enable_if_t<detector_t::materials::template is_defined<material_map<
                         rectangle2D<>, typename detector_t::scalar_type>>(),
                     void>> {
    using type = material_map<rectangle2D<>, typename detector_t::scalar_type>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_rectangle2_map};
};

/// Check for a 3D cuboid volume material map
template <typename detector_t>
struct mat_map_info<
    io::material_id::cuboid3_map, detector_t,
    std::enable_if_t<detector_t::materials::template is_defined<material_map<
                         cuboid3D<>, typename detector_t::scalar_type>>(),
                     void>> {
    using type = material_map<cuboid3D<>, typename detector_t::scalar_type>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_cuboid3_map};
};

/// Check for a 2D cylindrical material map
template <typename detector_t>
struct mat_map_info<
    io::material_id::cylinder2_map, detector_t,
    std::enable_if_t<detector_t::materials::template is_defined<material_map<
                         cylinder2D<>, typename detector_t::scalar_type>>(),
                     void>> {
    using type = material_map<cylinder2D<>, typename detector_t::scalar_type>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_cylinder2_map};
};

/// Check for a 3D cylindrical volume material map
template <typename detector_t>
struct mat_map_info<
    io::material_id::cylinder3_map, detector_t,
    std::enable_if_t<detector_t::materials::template is_defined<material_map<
                         cylinder3D, typename detector_t::scalar_type>>(),
                     void>> {
    using type = material_map<cylinder3D, typename detector_t::scalar_type>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_cylinder3_map};
};
/// @}

}  // namespace detray::io::detail
