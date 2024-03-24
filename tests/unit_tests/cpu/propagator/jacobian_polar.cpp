/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/propagator/detail/jacobian_polar.hpp"

#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/ring2D.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/tracks.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

using algebra_t = test::algebra;
using point3 = test::point3;
using vector3 = test::vector3;
using transform3 = test::transform3;
using matrix_operator = test::matrix_operator;
template <std::size_t ROWS, std::size_t COLS>
using matrix_type = test::matrix<ROWS, COLS>;

const scalar isclose{1e-5f};

// This test polar2D coordinate
GTEST_TEST(detray_propagator, jacobian_polar2D) {

    using jac_engine = detail::jacobian_engine<polar2D<algebra_t>>;

    // Preparation work
    const vector3 z = {0.f, 0.f, 1.f};
    const vector3 x = {1.f, 0.f, 0.f};
    const point3 t = {2.f, 3.f, 4.f};
    const transform3 trf(t, z, x);
    const point3 global1 = {4.f, 7.f, 4.f};
    const vector3 mom = {1.f, 2.f, 3.f};
    const scalar time{0.1f};
    const scalar charge{-1.};

    const scalar r{2.f};
    mask<ring2D> rng{0u, 0.f, r};

    // Free track parameter
    const free_track_parameters<algebra_t> free_params(global1, time, mom,
                                                       charge);
    const auto free_vec1 = free_params.vector();

    const auto bound_vec =
        detail::free_to_bound_vector<polar2D<algebra_t>>(trf, free_vec1);
    const auto free_vec2 = detail::bound_to_free_vector(trf, rng, bound_vec);

    const matrix_operator m;

    // Check if the bound vector is correct
    ASSERT_NEAR(m.element(bound_vec, 0u, 0u), std::sqrt(20.f), isclose);
    ASSERT_NEAR(m.element(bound_vec, 1u, 0u), std::atan2(4.f, 2.f), isclose);
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
    const bound_matrix<algebra_t> J =
        jac_engine::free_to_bound_jacobian(trf, free_vec1) *
        jac_engine::bound_to_free_jacobian(trf, rng, bound_vec);

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
