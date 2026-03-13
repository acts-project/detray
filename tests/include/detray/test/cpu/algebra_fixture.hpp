/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/type_traits.hpp"

// Test include(s).
#include "detray/test/framework/fixture_base.hpp"
#include "detray/test/framework/types.hpp"

namespace detray::test {

/// Test case class, contains the basic test definitions and algebra types
class detray_algebra : public fixture_base<> {
    using base = fixture_base<>;

    public:
    /// Constructor
    using base::base;

    protected:
#if !DETRAY_ALGEBRA_VC_AOS
    template <std::size_t ROWS, std::size_t COLS>
    void matrix_test_ops_any_matrix() {
        // Test the set_product method.
        {
            dmatrix<algebra_t, ROWS, ROWS> m1;
            dmatrix<algebra_t, ROWS, COLS> m2;

            for (std::size_t i = 0; i < ROWS; ++i) {
                for (std::size_t j = 0; j < ROWS; ++j) {
                    detray::getter::element(m1, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * ROWS + j);
                }
            }

            for (std::size_t i = 0; i < ROWS; ++i) {
                for (std::size_t j = 0; j < COLS; ++j) {
                    detray::getter::element(m2, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * COLS + j);
                }
            }

            {
                dmatrix<algebra_t, ROWS, COLS> r1 = m1 * m2;
                dmatrix<algebra_t, ROWS, COLS> r2;
                detray::matrix::set_product(r2, m1, m2);

                for (std::size_t i = 0; i < ROWS; ++i) {
                    for (std::size_t j = 0; j < COLS; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }
        }

        // Test the set_product_right_transpose method.
        {
            dmatrix<algebra_t, ROWS, ROWS> m1;
            dmatrix<algebra_t, COLS, ROWS> m2;

            for (std::size_t i = 0; i < ROWS; ++i) {
                for (std::size_t j = 0; j < ROWS; ++j) {
                    detray::getter::element(m1, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * ROWS + j);
                }
            }

            for (std::size_t i = 0; i < COLS; ++i) {
                for (std::size_t j = 0; j < ROWS; ++j) {
                    detray::getter::element(m2, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * COLS + j);
                }
            }

            {
                dmatrix<algebra_t, ROWS, COLS> r1 =
                    m1 * detray::matrix::transpose(m2);
                dmatrix<algebra_t, ROWS, COLS> r2;
                detray::matrix::set_product_right_transpose(r2, m1, m2);

                for (std::size_t i = 0; i < ROWS; ++i) {
                    for (std::size_t j = 0; j < COLS; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }
        }

        // Test the set_product_left_transpose method.
        {
            dmatrix<algebra_t, ROWS, ROWS> m1;
            dmatrix<algebra_t, ROWS, COLS> m2;

            for (std::size_t i = 0; i < ROWS; ++i) {
                for (std::size_t j = 0; j < ROWS; ++j) {
                    detray::getter::element(m1, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * ROWS + j);
                }
            }

            for (std::size_t i = 0; i < ROWS; ++i) {
                for (std::size_t j = 0; j < COLS; ++j) {
                    detray::getter::element(m2, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * COLS + j);
                }
            }

            {
                dmatrix<algebra_t, ROWS, COLS> r1 =
                    detray::matrix::transpose(m1) * m2;
                dmatrix<algebra_t, ROWS, COLS> r2;
                detray::matrix::set_product_left_transpose(r2, m1, m2);

                for (std::size_t i = 0; i < ROWS; ++i) {
                    for (std::size_t j = 0; j < COLS; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }

            // Test the transposable_product method.
            {
                // Both untransposed
                {
                    dmatrix<algebra_t, ROWS, COLS> r1 = m1 * m2;
                    dmatrix<algebra_t, ROWS, COLS> r2 =
                        detray::matrix::transposed_product<false, false>(m1,
                                                                         m2);

                    for (std::size_t i = 0; i < ROWS; ++i) {
                        for (std::size_t j = 0; j < COLS; ++j) {
                            ASSERT_NEAR(detray::getter::element(r1, i, j),
                                        detray::getter::element(r2, i, j),
                                        this->epsilon());
                        }
                    }
                }
            }
        }
    }

    template <std::size_t N>
    void matrix_test_ops_square_matrix() {
        {
            dmatrix<algebra_t, N, N> m1;
            dmatrix<algebra_t, N, N> m2;

            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < N; ++j) {
                    detray::getter::element(m1, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * N + j);
                    detray::getter::element(m2, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            -1 * (i * N + j) + 42);
                }
            }

            // Test the set_product method.
            {
                dmatrix<algebra_t, N, N> r1 = m1 * m2;
                dmatrix<algebra_t, N, N> r2;
                detray::matrix::set_product(r2, m1, m2);

                for (std::size_t i = 0; i < N; ++i) {
                    for (std::size_t j = 0; j < N; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }

            // Test the set_product_right_transpose method.
            {
                dmatrix<algebra_t, N, N> r1 =
                    m1 * detray::matrix::transpose(m2);
                dmatrix<algebra_t, N, N> r2;
                detray::matrix::set_product_right_transpose(r2, m1, m2);

                for (std::size_t i = 0; i < N; ++i) {
                    for (std::size_t j = 0; j < N; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }

            // Test the set_product_left_transpose method.
            {
                dmatrix<algebra_t, N, N> r1 =
                    detray::matrix::transpose(m1) * m2;
                dmatrix<algebra_t, N, N> r2;
                detray::matrix::set_product_left_transpose(r2, m1, m2);

                for (std::size_t i = 0; i < N; ++i) {
                    for (std::size_t j = 0; j < N; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }

            // Test the set_inplace_product_right method.
            {
                dmatrix<algebra_t, N, N> r1 = m1 * m2;
                dmatrix<algebra_t, N, N> r2 = m1;
                detray::matrix::set_inplace_product_right(r2, m2);

                for (std::size_t i = 0; i < N; ++i) {
                    for (std::size_t j = 0; j < N; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }

            // Test the set_inplace_product_left method.
            {
                dmatrix<algebra_t, N, N> r1 = m1 * m2;
                dmatrix<algebra_t, N, N> r2 = m2;
                detray::matrix::set_inplace_product_left(r2, m1);

                for (std::size_t i = 0; i < N; ++i) {
                    for (std::size_t j = 0; j < N; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }

            // Test the set_inplace_product_right_transpose method.
            {
                dmatrix<algebra_t, N, N> r1 =
                    m1 * detray::matrix::transpose(m2);
                dmatrix<algebra_t, N, N> r2 = m1;
                detray::matrix::set_inplace_product_right_transpose(r2, m2);

                for (std::size_t i = 0; i < N; ++i) {
                    for (std::size_t j = 0; j < N; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }

            // Test the set_inplace_product_left_transpose method.
            {
                dmatrix<algebra_t, N, N> r1 =
                    detray::matrix::transpose(m1) * m2;
                dmatrix<algebra_t, N, N> r2 = m2;
                detray::matrix::set_inplace_product_left_transpose(r2, m1);

                for (std::size_t i = 0; i < N; ++i) {
                    for (std::size_t j = 0; j < N; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }

            // Test the transposable_product method.
            {
                // Only left transposed
                {
                    dmatrix<algebra_t, N, N> r1 =
                        detray::matrix::transpose(m1) * m2;
                    dmatrix<algebra_t, N, N> r2 =
                        detray::matrix::transposed_product<true, false>(m1, m2);

                    for (std::size_t i = 0; i < N; ++i) {
                        for (std::size_t j = 0; j < N; ++j) {
                            ASSERT_NEAR(detray::getter::element(r1, i, j),
                                        detray::getter::element(r2, i, j),
                                        this->epsilon());
                        }
                    }
                }

                // Only right transposed
                {
                    dmatrix<algebra_t, N, N> r1 =
                        m1 * detray::matrix::transpose(m2);
                    dmatrix<algebra_t, N, N> r2 =
                        detray::matrix::transposed_product<false, true>(m1, m2);

                    for (std::size_t i = 0; i < N; ++i) {
                        for (std::size_t j = 0; j < N; ++j) {
                            ASSERT_NEAR(detray::getter::element(r1, i, j),
                                        detray::getter::element(r2, i, j),
                                        this->epsilon());
                        }
                    }
                }

                // Both transposed
                {
                    dmatrix<algebra_t, N, N> r1 =
                        detray::matrix::transpose(m1) *
                        detray::matrix::transpose(m2);
                    dmatrix<algebra_t, N, N> r2 =
                        detray::matrix::transposed_product<true, true>(m1, m2);

                    for (std::size_t i = 0; i < N; ++i) {
                        for (std::size_t j = 0; j < N; ++j) {
                            ASSERT_NEAR(detray::getter::element(r1, i, j),
                                        detray::getter::element(r2, i, j),
                                        this->epsilon());
                        }
                    }
                }
            }
        }

        this->template matrix_test_ops_any_matrix<N, N>();
    }

    template <std::size_t M, std::size_t N, std::size_t O>
    void matrix_test_ops_inhomogeneous_multipliable_matrices() {
        // Test NxM and MxO matrix multiplication
        {
            dmatrix<algebra_t, N, M> m1;
            dmatrix<algebra_t, M, O> m2;

            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < M; ++j) {
                    detray::getter::element(m1, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * N + j);
                }
            }

            for (std::size_t i = 0; i < M; ++i) {
                for (std::size_t j = 0; j < O; ++j) {
                    detray::getter::element(m2, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * M + j);
                }
            }

            {
                dmatrix<algebra_t, N, O> r1 = m1 * m2;
                dmatrix<algebra_t, N, O> r2 =
                    detray::matrix::transposed_product<false, false>(m1, m2);

                for (std::size_t i = 0; i < N; ++i) {
                    for (std::size_t j = 0; j < O; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }
        }

        // Test NxM and (OxM)^T matrix multiplication
        {
            dmatrix<algebra_t, N, M> m1;
            dmatrix<algebra_t, O, M> m2;

            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < M; ++j) {
                    detray::getter::element(m1, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * N + j);
                }
            }

            for (std::size_t i = 0; i < O; ++i) {
                for (std::size_t j = 0; j < M; ++j) {
                    detray::getter::element(m2, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * O + j);
                }
            }

            {
                dmatrix<algebra_t, N, O> r1 =
                    m1 * detray::matrix::transpose(m2);
                dmatrix<algebra_t, N, O> r2 =
                    detray::matrix::transposed_product<false, true>(m1, m2);

                for (std::size_t i = 0; i < N; ++i) {
                    for (std::size_t j = 0; j < O; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }
        }

        // Test (MxN)^T and MxO matrix multiplication
        {
            dmatrix<algebra_t, M, N> m1;
            dmatrix<algebra_t, M, O> m2;

            for (std::size_t i = 0; i < M; ++i) {
                for (std::size_t j = 0; j < N; ++j) {
                    detray::getter::element(m1, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * M + j);
                }
            }

            for (std::size_t i = 0; i < M; ++i) {
                for (std::size_t j = 0; j < O; ++j) {
                    detray::getter::element(m2, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * M + j);
                }
            }

            {
                dmatrix<algebra_t, N, O> r1 =
                    detray::matrix::transpose(m1) * m2;
                dmatrix<algebra_t, N, O> r2 =
                    detray::matrix::transposed_product<true, false>(m1, m2);

                for (std::size_t i = 0; i < N; ++i) {
                    for (std::size_t j = 0; j < O; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }
        }

        // Test (MxN)^T and (OxM)^T matrix multiplication
        {
            dmatrix<algebra_t, M, N> m1;
            dmatrix<algebra_t, O, M> m2;

            for (std::size_t i = 0; i < M; ++i) {
                for (std::size_t j = 0; j < N; ++j) {
                    detray::getter::element(m1, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * M + j);
                }
            }

            for (std::size_t i = 0; i < O; ++i) {
                for (std::size_t j = 0; j < M; ++j) {
                    detray::getter::element(m2, i, j) =
                        static_cast<detray::traits::scalar_t<decltype(m1)>>(
                            i * O + j);
                }
            }

            {
                dmatrix<algebra_t, N, O> r1 = detray::matrix::transpose(m1) *
                                              detray::matrix::transpose(m2);
                dmatrix<algebra_t, N, O> r2 =
                    detray::matrix::transposed_product<true, true>(m1, m2);

                for (std::size_t i = 0; i < N; ++i) {
                    for (std::size_t j = 0; j < O; ++j) {
                        ASSERT_NEAR(detray::getter::element(r1, i, j),
                                    detray::getter::element(r2, i, j),
                                    this->epsilon());
                    }
                }
            }
        }
    }
#endif
};

}  // namespace detray::test
