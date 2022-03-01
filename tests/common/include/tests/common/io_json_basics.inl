/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/io/json_io.hpp"

/// @note __plugin has to be defined with a preprocessor command

TEST(IO_JSON, write_read_transform) {

    detray::scalar x = 10.;
    detray::scalar y = 20.;
    detray::scalar z = 30.;

    // Let's make a matrix for
    detray::matrix<detray::scalar, 4, 4> ma;
    algebra::getter::element(ma, 0, 0) = 0.36;
    algebra::getter::element(ma, 0, 1) = 0.48;
    algebra::getter::element(ma, 0, 2) = -0.80;
    algebra::getter::element(ma, 0, 3) = x;
    algebra::getter::element(ma, 1, 0) = -0.80;
    algebra::getter::element(ma, 1, 1) = 0.60;
    algebra::getter::element(ma, 1, 2) = 0.;
    algebra::getter::element(ma, 1, 3) = y;
    algebra::getter::element(ma, 2, 0) = 0.48;
    algebra::getter::element(ma, 2, 1) = 0.64;
    algebra::getter::element(ma, 2, 2) = 0.60;
    algebra::getter::element(ma, 2, 3) = z;
    algebra::getter::element(ma, 3, 0) = 0.;
    algebra::getter::element(ma, 3, 1) = 0.;
    algebra::getter::element(ma, 3, 2) = 0.;
    algebra::getter::element(ma, 3, 2) = 1.;

    __plugin::transform3<detray::scalar> t3_out(ma);

    // Create the json command for the
    nlohmann::json tj = t3_out;

    nlohmann::json j;
    j["transform"] = tj;

    auto t3_in = j["transform"].get<__plugin::transform3<detray::scalar>>();

    auto translation_in = t3_in.translation();

    constexpr detray::scalar epsilon =
        std::numeric_limits<detray::scalar>::epsilon();

    ASSERT_NEAR(translation_in[0], x, epsilon);
    ASSERT_NEAR(translation_in[1], y, epsilon);
    ASSERT_NEAR(translation_in[2], z, epsilon);
}
