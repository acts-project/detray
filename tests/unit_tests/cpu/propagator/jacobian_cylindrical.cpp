/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/propagator/detail/jacobian_cylindrical.hpp"

#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/cylinder2D.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/tracks.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <limits>

using namespace detray;

using point3 = test::point3;
using vector3 = test::vector3;
using transform3 = test::transform3;
using matrix_operator = typename transform3::matrix_actor;
using size_type = typename matrix_operator::size_ty;
template <size_type ROWS, size_type COLS>
using matrix_type = typename matrix_operator::template matrix_type<ROWS, COLS>;

constexpr scalar isclose{1e-5f};

// This test cylindrical2D coordinate
GTEST_TEST(detray_propagator, jacobian_cylindrical2D) {

    using jac_engine = detail::jacobian_engine<cylindrical2D<transform3>>;

    // Preparation work
    const vector3 z = {0.f, 0.f, 1.f};
    const vector3 x = {1.f, 0.f, 0.f};
    const point3 t = {2.f, 3.f, 4.f};
    const transform3 trf(t, z, x);
    // Global position on surface
    const point3 global1 = {3.4142136f, 4.4142136f, 9.f};
    const vector3 mom = {1.f, 2.f, 3.f};
    const scalar time{0.1f};
    const scalar charge{-1.f};

    const scalar r{2.f};
    const scalar hz{detail::invalid_value<scalar>()};
    mask<cylinder2D> cyl{0u, r, -hz, hz};

    // Free track parameter
    const free_track_parameters<transform3> free_params(global1, time, mom,
                                                        charge);
    const auto free_vec1 = free_params.vector();

    const auto bound_vec =
        detail::free_to_bound_vector<cylindrical2D<transform3>>(trf, free_vec1);
    const auto free_vec2 = detail::bound_to_free_vector(trf, cyl, bound_vec);

    const matrix_operator m;

    // Check if the bound vector is correct
    ASSERT_NEAR(m.element(bound_vec, 0u, 0u), r * constant<scalar>::pi_4,
                isclose);
    ASSERT_NEAR(m.element(bound_vec, 1u, 0u), 5.f, isclose);
    ASSERT_NEAR(m.element(bound_vec, 2u, 0u), 1.1071487f,
                isclose);  // atan(2)
    ASSERT_NEAR(m.element(bound_vec, 3u, 0u), 0.64052231f,
                isclose);  // atan(sqrt(5)/3)
    ASSERT_NEAR(m.element(bound_vec, 4u, 0u), -1.f / 3.7416574f, isclose);
    ASSERT_NEAR(m.element(bound_vec, 5u, 0u), 0.1f, isclose);

    // Check if the same free vector is obtained
    for (unsigned int i = 0u; i < 8u; i++) {
        ASSERT_NEAR(m.element(free_vec1, i, 0u), m.element(free_vec2, i, 0u),
                    isclose);
    }

    // Test Jacobian transformation
    const bound_matrix<transform3> J =
        jac_engine::free_to_bound_jacobian(trf, free_vec1) *
        jac_engine::bound_to_free_jacobian(trf, cyl, bound_vec);

    for (unsigned int i = 0u; i < 6u; i++) {
        for (unsigned int j = 0u; j < 6u; j++) {
            if (i == j) {
                EXPECT_NEAR(m.element(J, i, j), 1.f, isclose);
            } else {
                EXPECT_NEAR(m.element(J, i, j), 0.f, isclose);
            }
        }
    }
}
