/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#if DETRAY_ALGEBRA_ARRAY
#include "detray/plugins/algebra/array_definitions.hpp"
#elif DETRAY_ALGEBRA_EIGEN
#include "detray/plugins/algebra/eigen_definitions.hpp"
#elif DETRAY_ALGEBRA_SMATRIX
#include "detray/plugins/algebra/smatrix_definitions.hpp"
#elif DETRAY_ALGEBRA_VC
#include "detray/plugins/algebra/vc_array_definitions.hpp"
#elif DETRAY_ALGEBRA_VC_SOA
#include "detray/plugins/algebra/vc_soa_definitions.hpp"
#else
#error "No algebra plugin selected! Please link to one of the algebra plugins."
#endif

// System include(s)
#include <type_traits>

namespace detray {

namespace detail {
/// The detray scalar types (can be SIMD)
/// @{
template <typename T, typename = void>
struct get_scalar {};

template <typename T>
struct get_scalar<T, std::enable_if_t<std::is_arithmetic_v<T>, void>> {
    using scalar = T;
};

template <typename T>
struct get_scalar<
    T, std::enable_if_t<!std::is_same_v<typename T::scalar, void>, void>> {
    using scalar = typename T::scalar;
};
/// @}

/// The detray algebra types (can be SIMD)
/// @{
template <typename T, typename = void>
struct get_algebra {};

template <typename T>
struct get_algebra<T, std::enable_if_t<std::is_arithmetic_v<T>, void>> {
    // vectors and transforms defined in 4D homogeneous coordinates
    using point3D = std::array<T, 4>;
    using vector3D = std::array<T, 4>;
    using transform3D = std::array<T, 16>;
};

template <typename T>
struct get_algebra<
    T, std::enable_if_t<!std::is_same_v<typename T::point3D, void>, void>> {
    using point2D = typename T::point2D;
    using point3D = typename T::point3D;
    using vector3D = typename T::vector3D;
    using transform3D = typename T::transform3D;
};
/// @}

/// The detray matrix types
/// @{
template <typename T, typename = void>
struct get_matrix {};

template <typename T>
struct get_matrix<
    T, std::enable_if_t<!std::is_same_v<typename T::matrix_operator, void>,
                        void>> {
    using matrix_operator = typename T::matrix_operator;
    using size_type = typename matrix_operator::size_ty;

    template <std::size_t ROWS, std::size_t COLS>
    using matrix = typename matrix_operator::template matrix_type<
        static_cast<size_type>(ROWS), static_cast<size_type>(COLS)>;
};
/// @}

}  // namespace detail

template <template <typename> class A, typename T>
using dsimd = typename A<float>::template simd<T>;

template <typename A = detray::scalar>
using dscalar = typename detail::get_scalar<A>::scalar;

template <typename A>
using dpoint2D = typename detail::get_algebra<A>::point2D;

template <typename A>
using dpoint3D = typename detail::get_algebra<A>::point3D;

template <typename A>
using dvector3D = typename detail::get_algebra<A>::vector3D;

template <typename A>
using dtransform3D = typename detail::get_algebra<A>::transform3D;

template <typename A>
using dmatrix_operator = typename detail::get_matrix<A>::matrix_operator;

template <typename A>
using dsize_type = typename detail::get_matrix<A>::size_type;

template <typename A, std::size_t R, std::size_t C>
using dmatrix = typename detail::get_matrix<A>::template matrix<R, C>;

namespace detail {

/// Check if an algebra has soa layout
/// @{
template <typename A, typename = void>
struct is_soa : public std::false_type {};

template <typename A>
struct is_soa<A, std::enable_if_t<!std::is_arithmetic_v<dscalar<A>>, void>>
    : public std::true_type {};

template <typename A>
inline constexpr bool is_soa_v = is_soa<A>::value;
/// @}

}  // namespace detail

}  // namespace detray
