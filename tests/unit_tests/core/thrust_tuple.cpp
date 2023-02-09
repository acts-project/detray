/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Thrust include(s).
#include <thrust/tuple.h>

// GoogleTest include(s).
#include <gtest/gtest.h>

namespace {

struct TestStruct {

    int m_int = 0;
    float m_float = 0.f;
    double m_double = 0.;

};  // struct TestStruct

}  // namespace

TEST(thrust_tuple, constructor) {

    auto t1 = thrust::make_tuple(12, 2.f, 3.);
    (void)t1;
    auto t2 = thrust::make_tuple(3.14f, TestStruct{3, 4.f, 5.});
    (void)t2;
}

TEST(thrust_tuple, access) {

    auto t1 = thrust::make_tuple(3.14f, TestStruct{3, 4.f, 5.});
    EXPECT_FLOAT_EQ(thrust::get<0>(t1), 3.14f);
    EXPECT_EQ(thrust::get<1>(t1).m_int, 3);
    EXPECT_FLOAT_EQ(thrust::get<1>(t1).m_float, 4.f);
    EXPECT_DOUBLE_EQ(thrust::get<1>(t1).m_double, 5.);
}
