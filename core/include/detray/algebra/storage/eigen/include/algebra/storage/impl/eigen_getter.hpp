/** Algebra plugins, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/concepts.hpp"
#include "algebra/math/common.hpp"
#include "algebra/qualifiers.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 20012
#endif  // __NVCC_DIAG_PRAGMA_SUPPORT__
#include <Eigen/Core>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic pop
#endif  // __NVCC_DIAG_PRAGMA_SUPPORT__

// System include(s).
#include <type_traits>

namespace algebra::eigen::storage {

/// Functor used to access elements of Eigen matrices
struct element_getter {
  /// Get non-const access to a matrix element
  template <typename derived_type, concepts::index size_type_1,
            concepts::index size_type_2>
  requires std::is_base_of_v<
      Eigen::DenseCoeffsBase<derived_type, Eigen::WriteAccessors>,
      Eigen::MatrixBase<derived_type> >
      ALGEBRA_HOST_DEVICE constexpr auto &operator()(
          Eigen::MatrixBase<derived_type> &m, size_type_1 row,
          size_type_2 col) const {

    return m(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
  }
  /// Get const access to a matrix element
  template <typename derived_type, concepts::index size_type_1,
            concepts::index size_type_2>
  ALGEBRA_HOST_DEVICE constexpr auto operator()(
      const Eigen::MatrixBase<derived_type> &m, size_type_1 row,
      size_type_2 col) const {

    return m(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
  }
  /// Get non-const access to a matrix element
  template <typename derived_type, concepts::index size_type>
  requires std::is_base_of_v<
      Eigen::DenseCoeffsBase<derived_type, Eigen::WriteAccessors>,
      Eigen::MatrixBase<derived_type> >
      ALGEBRA_HOST_DEVICE constexpr auto &operator()(
          Eigen::MatrixBase<derived_type> &m, size_type row) const {

    return m(static_cast<Eigen::Index>(row));
  }
  /// Get const access to a matrix element
  template <typename derived_type, concepts::index size_type>
  ALGEBRA_HOST_DEVICE constexpr auto operator()(
      const Eigen::MatrixBase<derived_type> &m, size_type row) const {

    return m(static_cast<Eigen::Index>(row));
  }
};  // struct element_getter

/// Function extracting an element from a matrix (const)
template <typename derived_type>
ALGEBRA_HOST_DEVICE constexpr decltype(auto) element(
    const Eigen::MatrixBase<derived_type> &m, std::size_t row,
    std::size_t col) {

  return element_getter()(m, static_cast<Eigen::Index>(row),
                          static_cast<Eigen::Index>(col));
}

/// Function extracting an element from a matrix (non-const)
template <typename derived_type>
requires std::is_base_of_v<
    Eigen::DenseCoeffsBase<derived_type, Eigen::WriteAccessors>,
    Eigen::MatrixBase<derived_type> >
    ALGEBRA_HOST_DEVICE constexpr decltype(auto) element(
        Eigen::MatrixBase<derived_type> &m, std::size_t row, std::size_t col) {

  return element_getter()(m, static_cast<Eigen::Index>(row),
                          static_cast<Eigen::Index>(col));
}

/// Function extracting an element from a matrix (const)
template <typename derived_type>
ALGEBRA_HOST_DEVICE constexpr decltype(auto) element(
    const Eigen::MatrixBase<derived_type> &m, std::size_t row) {

  return element_getter()(m, static_cast<Eigen::Index>(row));
}

/// Function extracting an element from a matrix (non-const)
template <typename derived_type>
requires std::is_base_of_v<
    Eigen::DenseCoeffsBase<derived_type, Eigen::WriteAccessors>,
    Eigen::MatrixBase<derived_type> >
    ALGEBRA_HOST_DEVICE constexpr decltype(auto) element(
        Eigen::MatrixBase<derived_type> &m, std::size_t row) {

  return element_getter()(m, static_cast<Eigen::Index>(row));
}

/// Functor used to extract a block from Eigen matrices
struct block_getter {

  template <int kROWS, int kCOLS, typename derived_type,
            concepts::index size_type_1, concepts::index size_type_2>
  ALGEBRA_HOST_DEVICE decltype(auto) operator()(
      const Eigen::MatrixBase<derived_type> &m, size_type_1 row,
      size_type_2 col) const {

    return m.template block<kROWS, kCOLS>(row, col);
  }

  template <int kROWS, int kCOLS, typename derived_type,
            concepts::index size_type_1, concepts::index size_type_2>
  ALGEBRA_HOST_DEVICE decltype(auto) operator()(
      Eigen::MatrixBase<derived_type> &m, size_type_1 row,
      size_type_2 col) const {

    return m.template block<kROWS, kCOLS>(row, col);
  }

  template <int SIZE, typename derived_type, concepts::index size_type_1,
            concepts::index size_type_2>
  ALGEBRA_HOST_DEVICE decltype(auto) vector(Eigen::MatrixBase<derived_type> &m,
                                            size_type_1 row,
                                            size_type_2 col) const {

    return m.template block<SIZE, 1>(row, col);
  }

  template <int SIZE, typename derived_type, concepts::index size_type_1,
            concepts::index size_type_2>
  ALGEBRA_HOST_DEVICE decltype(auto) vector(
      const Eigen::MatrixBase<derived_type> &m, size_type_1 row,
      size_type_2 col) const {

    return m.template block<SIZE, 1>(row, col);
  }

};  // struct block_getter

/// Operator getting a block of a const matrix
template <int ROWS, int COLS, class derived_type>
ALGEBRA_HOST_DEVICE decltype(auto) block(
    const Eigen::MatrixBase<derived_type> &m, std::size_t row,
    std::size_t col) {
  return block_getter{}.template operator()<ROWS, COLS>(
      m, static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
}

/// Operator getting a block of a const matrix
template <int ROWS, int COLS, class derived_type>
ALGEBRA_HOST_DEVICE decltype(auto) block(Eigen::MatrixBase<derived_type> &m,
                                         std::size_t row, std::size_t col) {
  return block_getter{}.template operator()<ROWS, COLS>(
      m, static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
}

/// Function extracting a slice from the matrix
template <int SIZE, typename derived_type>
ALGEBRA_HOST_DEVICE constexpr decltype(auto) vector(
    const Eigen::MatrixBase<derived_type> &m, std::size_t row,
    std::size_t col) {

  return block_getter{}.template vector<SIZE>(m, static_cast<Eigen::Index>(row),
                                              static_cast<Eigen::Index>(col));
}

/// Operator setting a block
template <typename derived_type1, typename derived_type2>
ALGEBRA_HOST_DEVICE void set_block(Eigen::MatrixBase<derived_type1> &m,
                                   const Eigen::MatrixBase<derived_type2> &b,
                                   std::size_t row, std::size_t col) {
  using block_t = Eigen::MatrixBase<derived_type2>;
  constexpr auto R{block_t::RowsAtCompileTime};
  constexpr auto C{block_t::ColsAtCompileTime};
  m.template block<R, C>(static_cast<Eigen::Index>(row),
                         static_cast<Eigen::Index>(col)) = b;
}

}  // namespace algebra::eigen::storage
