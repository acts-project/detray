/** Detray library, part of the ACTS project
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/masks/masks.hpp"
#include "detray/tracks/bound_track_parameters.hpp"

// Google test include(s).
#include <gtest/gtest.h>

using namespace detray;
using transform3_t = __plugin::transform3<detray::scalar>;

/// This tests the basic functionality of projection matrix
TEST(projection_matrix, rectangle) {

    bound_track_parameters<transform3_t> bound_params;

    // 2D, normal ordering
    mask<rectangle2D<>> rec_0{0u, 0.f, 0.f};

    // Check the projection matrix
    // [ 1 0 0 0 0 0 ]
    // [ 0 1 0 0 0 0 ]
    const auto proj_0 = rec_0.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj_0, i, j), 1.f);
            } else {
                ASSERT_EQ(getter::element(proj_0, i, j), 0.f);
            }
        }
    }

    // 2D, reverse ordering
    mask<rectangle2D<plane_intersector, 2u, false>> rec_1{0u, 0.f, 0.f};

    // Check the projection matrix
    // [ 0 1 0 0 0 0 ]
    // [ 1 0 0 0 0 0 ]
    const auto proj_1 = rec_1.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
                ASSERT_EQ(getter::element(proj_1, i, j), 1.f);
            } else {
                ASSERT_EQ(getter::element(proj_1, i, j), 0.f);
            }
        }
    }

    // 1D, normal ordering
    mask<rectangle2D<plane_intersector, 1u>> rec_2{0u, 0.f, 0.f};

    // Check the projection matrix
    // [ 1 0 0 0 0 0 ]
    const auto proj_2 = rec_2.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < 1u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj_2, i, j), 1.f);
            } else {
                ASSERT_EQ(getter::element(proj_2, i, j), 0.f);
            }
        }
    }

    // 1D, reverse ordering
    mask<rectangle2D<plane_intersector, 1u, false>> rec_3{0u, 0.f, 0.f};

    // Check the projection matrix
    // [ 0 1 0 0 0 0 ]
    const auto proj_3 = rec_3.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < 1u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == 0 && j == 1) {
                ASSERT_EQ(getter::element(proj_3, i, j), 1.f);
            } else {
                ASSERT_EQ(getter::element(proj_3, i, j), 0.f);
            }
        }
    }
}

TEST(projection_matrix, annulus) {

    bound_track_parameters<transform3_t> bound_params;

    // 1D reverse ordering
    mask<annulus2D<>> ann_0{0u, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};

    // Check the projection matrix
    // [ 0 1 0 0 0 0 ]
    const auto proj_0 = ann_0.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < 1u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == 0 && j == 1) {
                ASSERT_EQ(getter::element(proj_0, i, j), 1.f);
            } else {
                ASSERT_EQ(getter::element(proj_0, i, j), 0.f);
            }
        }
    }

    // 1D normal ordering
    mask<annulus2D<plane_intersector, 1u, true>> ann_1{0u,  0.f, 0.f, 0.f,
                                                       0.f, 0.f, 0.f, 0.f};

    // Check the projection matrix
    // [ 1 0 0 0 0 0 ]
    const auto proj_1 = ann_1.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < 1u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj_1, i, j), 1.f);
            } else {
                ASSERT_EQ(getter::element(proj_1, i, j), 0.f);
            }
        }
    }

    // 2D reverse ordering
    mask<annulus2D<plane_intersector, 2u>> ann_2{0u,  0.f, 0.f, 0.f,
                                                 0.f, 0.f, 0.f, 0.f};

    // Check the projection matrix
    // [ 0 1 0 0 0 0 ]
    // [ 1 0 0 0 0 0 ]
    const auto proj_2 = ann_2.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
                ASSERT_EQ(getter::element(proj_2, i, j), 1.f);
            } else {
                ASSERT_EQ(getter::element(proj_2, i, j), 0.f);
            }
        }
    }

    // 2D normal ordering
    mask<annulus2D<plane_intersector, 2u, true>> ann_3{0u,  0.f, 0.f, 0.f,
                                                       0.f, 0.f, 0.f, 0.f};

    // Check the projection matrix
    // [ 1 0 0 0 0 0 ]
    // [ 0 1 0 0 0 0 ]
    const auto proj_3 = ann_3.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj_3, i, j), 1.f);
            } else {
                ASSERT_EQ(getter::element(proj_3, i, j), 0.f);
            }
        }
    }
}

TEST(projection_matrix, line) {

    bound_track_parameters<transform3_t> bound_params;
    auto& bound_vec = bound_params.vector();
    getter::element(bound_vec, e_bound_loc0, 0u) = -1.f;

    // 1D normal ordering
    mask<line<>> ln_0{0u, 0.f, 0.f};

    // Check the projection matrix
    // [ -1 0 0 0 0 0 ]
    const auto proj_0 = ln_0.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < 1u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj_0, i, j), -1.f);
            } else {
                ASSERT_EQ(getter::element(proj_0, i, j), 0.f);
            }
        }
    }

    // 1D reverse ordering
    mask<line<false, line_intersector, 1u, false>> ln_1{0u, 0.f, 0.f};

    // Check the projection matrix
    // [ 0 1 0 0 0 0 ]
    const auto proj_1 = ln_1.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < 1u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == 0 && j == 1) {
                ASSERT_EQ(getter::element(proj_1, i, j), 1.f);
            } else {
                ASSERT_EQ(getter::element(proj_1, i, j), 0.f);
            }
        }
    }

    // 2D normal ordering
    mask<line<false, line_intersector, 2u>> ln_2{0u, 0.f, 0.f};

    // Check the projection matrix
    // [-1 0 0 0 0 0 ]
    // [ 0 1 0 0 0 0 ]
    const auto proj_2 = ln_2.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == 0 && j == 0) {
                ASSERT_EQ(getter::element(proj_2, i, j), -1.f);
            } else if (i == 1 && j == 1) {
                ASSERT_EQ(getter::element(proj_2, i, j), 1.f);
            } else {
                ASSERT_EQ(getter::element(proj_2, i, j), 0.f);
            }
        }
    }

    // 2D reverse ordering
    mask<line<false, line_intersector, 2u, false>> ln_3{0u, 0.f, 0.f};

    // Check the projection matrix
    // [ 0 1 0 0 0 0 ]
    // [-1 0 0 0 0 0 ]
    const auto proj_3 = ln_3.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == 0 && j == 1) {
                ASSERT_EQ(getter::element(proj_3, i, j), 1.f);
            } else if (i == 1 && j == 0) {
                ASSERT_EQ(getter::element(proj_3, i, j), -1.f);
            } else {
                ASSERT_EQ(getter::element(proj_3, i, j), 0.f);
            }
        }
    }
}
