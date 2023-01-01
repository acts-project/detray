/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/utils/tuple_helpers.hpp"

// System include(s)
#include <cassert>
#include <string>
#include <tuple>
#include <type_traits>

// Thrust include(s).
#include <thrust/tuple.h>

// GoogleTest include(s).
#include <gtest/gtest.h>

TEST(utils, tuple_helpers) {

    using namespace detray;

    // std::tuple test
    auto s_tuple = detail::make_tuple<std::tuple>(
        2.0f, -3L, std::string("std::tuple"), 4UL);
    static_assert(std::is_same_v<std::tuple<float, long, std::string, unsigned long>,
                                 decltype(s_tuple)>,
                  "detail::make_tuple failed for std::tuple");

    static_assert(std::is_same_v<detail::tuple_element_t<2, decltype(s_tuple)>,
                                 std::string>,
                  "detail::tuple_element retrieval failed for std::tuple");

    const auto s_tuple_size = detail::tuple_size_v<decltype(s_tuple)>;
    EXPECT_EQ(s_tuple_size, 4UL);

    EXPECT_FLOAT_EQ(detail::get<0>(s_tuple), 2.0f);
    EXPECT_EQ(detail::get<1>(s_tuple), -3L);
    EXPECT_EQ(detail::get<2>(s_tuple), std::string("std::tuple"));
    EXPECT_EQ(detail::get<3>(s_tuple), 4UL);
    EXPECT_FLOAT_EQ(detail::get<float>(s_tuple), 2.0f);
    EXPECT_EQ(detail::get<long>(s_tuple), -3L);
    EXPECT_EQ(detail::get<std::string>(s_tuple), std::string("std::tuple"));
    EXPECT_EQ(detail::get<unsigned long>(s_tuple), 4UL);

    // thrust::tuple test
    auto t_tuple = detail::make_tuple<thrust::tuple>(
        1.0f, 2UL, std::string("thrust::tuple"));
    static_assert(std::is_same_v<thrust::tuple<float, unsigned long, std::string>,
                                 decltype(t_tuple)>,
                  "detail::make_tuple failed for thrust::tuple");

    static_assert(
        std::is_same_v<detail::tuple_element_t<1, decltype(t_tuple)>, unsigned long>,
        "detail::tuple_element retrieval failed for thrust::tuple");

    const auto t_tuple_size = detail::tuple_size_v<decltype(t_tuple)>;
    EXPECT_EQ(t_tuple_size, 3UL);

    EXPECT_FLOAT_EQ(detail::get<0>(t_tuple), 1.0f);
    EXPECT_EQ(detail::get<1>(t_tuple), 2UL);
    EXPECT_EQ(detail::get<2>(t_tuple), std::string("thrust::tuple"));
    EXPECT_FLOAT_EQ(detail::get<float>(t_tuple), 1.0f);
    EXPECT_EQ(detail::get<unsigned long>(t_tuple), 2UL);
    EXPECT_EQ(detail::get<std::string>(t_tuple), std::string("thrust::tuple"));
}