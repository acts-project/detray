/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"
#include "algebra/storage/vc_aos.hpp"
#include "algebra/storage/vector.hpp"

// Vc include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Vc/Vc>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

// System include(s).
#include <concepts>
#include <utility>

namespace algebra::vc_aos::math {

/// This method retrieves phi from a vector @param v
template <algebra::concepts::vc_aos_vector vector_t>
ALGEBRA_HOST_DEVICE constexpr auto phi(const vector_t &v) {
  return algebra::math::atan2(v[1], v[0]);
}

/// This method retrieves the perpendicular magnitude of a vector @param v
template <algebra::concepts::vc_aos_vector vector_t>
ALGEBRA_HOST_DEVICE constexpr auto perp(const vector_t &v) {
  return algebra::math::sqrt(algebra::math::fma(v[0], v[0], v[1] * v[1]));
}

/// This method retrieves theta from a vector @param v
template <algebra::concepts::vc_aos_vector vector_t>
ALGEBRA_HOST_DEVICE constexpr auto theta(const vector_t &v) {
  return algebra::math::atan2(perp(v), v[2]);
}

/// Dot product between two input vectors
///
/// @tparam vector_t generic input vector type
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <algebra::concepts::vc_aos_vector vector_t1,
          algebra::concepts::vc_aos_vector vector_t2>
ALGEBRA_HOST_DEVICE constexpr auto dot(const vector_t1 &a, const vector_t2 &b) {
  return (a * b).sum();
}

/// This method retrieves the norm of a vector, no dimension restriction
///
/// @param v the input vector
template <algebra::concepts::vc_aos_vector vector_t>
ALGEBRA_HOST_DEVICE constexpr auto norm(const vector_t &v) {
  return algebra::math::sqrt(dot(v, v));
}

/// Get a normalized version of the input vector
///
/// @tparam vector_t generic input vector type
///
/// @param v the input vector
template <algebra::concepts::vc_aos_vector vector_t>
ALGEBRA_HOST_DEVICE constexpr auto normalize(const vector_t &v) {
  return v / norm(v);
}

/// This method retrieves the pseudo-rapidity from a vector or vector base with
/// rows >= 3
///
/// @param v the input vector
template <algebra::concepts::vc_aos_vector vector_t>
ALGEBRA_HOST_DEVICE constexpr auto eta(const vector_t &v) noexcept {
  return algebra::math::atanh(v[2] / norm(v));
}

/// Cross product between two input vectors - 3 Dim (single precision)
///
/// @tparam vector_t generic input vector type
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector representing the cross product
template <algebra::concepts::vc_aos_vector vector_t1,
          algebra::concepts::vc_aos_vector vector_t2>
requires(std::same_as<algebra::traits::value_t<vector_t1>, float>
             &&std::same_as<algebra::traits::value_t<vector_t2>, float>)
    ALGEBRA_HOST_DEVICE
    constexpr auto cross(const vector_t1 &a, const vector_t2 &b)
        -> decltype(a * b - b * a) {
  using T = algebra::traits::value_t<vector_t1>;
  using simd_array_t = Vc::SimdArray<T, 4>;

  // Mask to write on the last two elements: [0011]
  const typename simd_array_t::mask_type m{simd_array_t::IndexesFromZero() > 1};

  // a     = [a1, a2, a3, 0]
  // a_rot = [a2, a3, 0, a1]
  simd_array_t a_rot = static_cast<simd_array_t>(a).rotated(1);
  // a_rot = [a2, a3, a1, a2]
  a_rot(m) = a_rot.rotated(1);

  // Same for b
  simd_array_t b_rot = static_cast<simd_array_t>(b).rotated(1);
  b_rot(m) = b_rot.rotated(1);

  simd_array_t res = a_rot * b_rot.rotated(1) - a_rot.rotated(1) * b_rot;
  // Restore trailing zero
  res[3] = 0.f;

  return res;
}

/// Cross product between two input vectors - 3 Dim (double precision)
///
/// @tparam vector_t generic input vector type
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector representing the cross product
template <algebra::concepts::vc_aos_vector vector_t1,
          algebra::concepts::vc_aos_vector vector_t2>
requires(std::same_as<algebra::traits::value_t<vector_t1>, double> ||
         std::same_as<algebra::traits::value_t<vector_t2>, double>)
    ALGEBRA_HOST_DEVICE
    constexpr auto cross(const vector_t1 &a, const vector_t2 &b)
        -> decltype(a * b - b * a) {

  return {algebra::math::fma(a[1], b[2], -b[1] * a[2]),
          algebra::math::fma(a[2], b[0], -b[2] * a[0]),
          algebra::math::fma(a[0], b[1], -b[0] * a[1]), 0.f};
}

/// Elementwise sum
///
/// @tparam vector_t generic input vector type
///
/// @param v the vector whose elements should be summed
///
/// @return the sum of the elements
template <algebra::concepts::vc_aos_vector vector_t>
ALGEBRA_HOST_DEVICE constexpr auto sum(const vector_t &v) {
  return v.get().sum();
}

}  // namespace algebra::vc_aos::math
