/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// covfie core
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>

// GTest include(s)
#include <gtest/gtest.h>

GTEST_TEST(Covfie, ConstantField1D) {
    using field_t =
        covfie::field<covfie::backend::constant<covfie::vector::float1,
                                                covfie::vector::float1>>;

    field_t f(
        covfie::make_parameter_pack(field_t::backend_t::configuration_t{2.f}));
    field_t::view_t v(f);

    for (float x = -100.f; x <= 100.f; x += 1.f) {
        EXPECT_EQ(v.at(x)[0], 2.f);
    }
}

GTEST_TEST(Covfie, ConstantField2D) {
    using field_t =
        covfie::field<covfie::backend::constant<covfie::vector::float2,
                                                covfie::vector::float2>>;

    field_t f(covfie::make_parameter_pack(
        field_t::backend_t::configuration_t{2.f, 5.f}));
    field_t::view_t v(f);

    for (float x = -100.f; x <= 100.f; x += 1.f) {
        for (float y = -100.f; y <= 100.f; y += 1.f) {
            EXPECT_EQ(v.at(x, y)[0], 2.f);
            EXPECT_EQ(v.at(x, y)[1], 5.f);
        }
    }
}

GTEST_TEST(Covfie, ConstantField3D) {
    using field_t =
        covfie::field<covfie::backend::constant<covfie::vector::float3,
                                                covfie::vector::float3>>;

    field_t f(covfie::make_parameter_pack(
        field_t::backend_t::configuration_t{2.f, 5.f, -4.f}));
    field_t::view_t v(f);

    for (float x = -10.f; x <= 10.f; x += 1.f) {
        for (float y = -10.f; y <= 10.f; y += 1.f) {
            for (float z = -10.f; z <= 10.f; z += 1.f) {
                EXPECT_EQ(v.at(x, y, z)[0], 2.f);
                EXPECT_EQ(v.at(x, y, z)[1], 5.f);
                EXPECT_EQ(v.at(x, y, z)[2], -4.f);
            }
        }
    }
}
