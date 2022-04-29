/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// detray include(s)
#include "detray/definitions/units.hpp"
#include "detray/propagator/detail/covariance_engine.hpp"
#include "detray/propagator/track.hpp"

// google-test include(s)
#include <gtest/gtest.h>

using namespace detray;
using size_type = __plugin::size_type;
using vector3 = __plugin::vector3<scalar>;
using covariance_engine = detail::covariance_engine<scalar>;
using matrix_operator = covariance_engine::matrix_operator;
using jacobian_engine = typename covariance_engine::jacobian_engine;
using vector_engine = typename covariance_engine::vector_engine;
using transform3 = typename covariance_engine::transform3;

constexpr scalar epsilon = 1e-3;

TEST(covariance_engine, jacobian_coordinate) {

    // test surface
    const vector3 u{0, 1, 0};
    const vector3 w{1, 0, 0};
    const vector3 t{0, 0, 0};
    const transform3 trf(t, w, u);

    // Bound track parameter
    vector3 local{2, 3, 0};
    vector3 mom{1., 2.5, 3.};
    scalar time = 0.;
    scalar q = -1.;

    // bound vector
    typename bound_track_parameters::vector_type bound_vec;
    getter::element(bound_vec, e_bound_loc0, 0) = local[0];
    getter::element(bound_vec, e_bound_loc1, 0) = local[1];
    getter::element(bound_vec, e_bound_phi, 0) = getter::phi(mom);
    getter::element(bound_vec, e_bound_theta, 0) = getter::theta(mom);
    getter::element(bound_vec, e_bound_qoverp, 0) = q / getter::norm(mom);
    getter::element(bound_vec, e_bound_time, 0) = time;

    typename free_track_parameters::vector_type free_vec =
        vector_engine().bound_to_free_vector(trf, bound_vec);

    // Make sure that vector engine works fine
    typename bound_track_parameters::vector_type new_bound_vec =
        vector_engine().free_to_bound_vector(trf, free_vec);
    for (size_type i = 0; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec, i, 0),
                    matrix_operator().element(new_bound_vec, i, 0), epsilon);
    }

    // Get Jacobian for coordinate transform
    const auto bound_to_free_matrix =
        jacobian_engine().bound_to_free_coordinate(trf, bound_vec);

    const auto free_to_bound_matrix =
        jacobian_engine().free_to_bound_coordinate(trf, free_vec);

    const bound_matrix I_bb = free_to_bound_matrix * bound_to_free_matrix;

    for (size_type i = 0; i < e_bound_size; i++) {
        for (size_type j = 0; j < e_bound_size; j++) {
            if (i == j) {
                EXPECT_NEAR(matrix_operator().element(I_bb, i, j), 1, epsilon);
            } else {
                EXPECT_NEAR(matrix_operator().element(I_bb, i, j), 0, epsilon);
            }
        }
    }
}