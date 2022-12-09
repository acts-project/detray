/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/utils/tuple_helpers.hpp"

#include <iostream>

// Thrust include(s).
#include <thrust/tuple.h>

// GoogleTest include(s).
#include <gtest/gtest.h>

TEST(tuple_helpers, tuple_helpers) {

    using namespace detray;

    // std::tuple test
    auto s_tuple = detail::make_tuple<std::tuple>(1.0, 2, "std::tuple");

    const auto s_tuple_size = detail::tuple_size<decltype(s_tuple)>::value;

    EXPECT_EQ(s_tuple_size, 3);
    EXPECT_FLOAT_EQ(detail::get<0>(s_tuple), 1.0);
    EXPECT_EQ(detail::get<1>(s_tuple), 2);
    EXPECT_EQ(detail::get<2>(s_tuple), "std::tuple");

    // thrust::tuple test
    auto t_tuple = detail::make_tuple<thrust::tuple>(1.0, 2, "thrust::tuple");

    const auto t_tuple_size = detail::tuple_size<decltype(t_tuple)>::value;

    EXPECT_EQ(t_tuple_size, 3);
    EXPECT_FLOAT_EQ(detail::get<0>(t_tuple), 1.0);
    EXPECT_EQ(detail::get<1>(t_tuple), 2);
    EXPECT_EQ(detail::get<2>(t_tuple), "thrust::tuple");
}