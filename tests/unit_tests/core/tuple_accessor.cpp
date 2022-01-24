/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <iostream>

#include "detray/definitions/detail/accessor.hpp"

// Thrust include(s).
#include <thrust/tuple.h>

// GoogleTest include(s).
#include <gtest/gtest.h>

namespace {

struct TestStruct {

    int m_int = 0;
    float m_float = 0.;
    double m_double = 0.;

};  // struct TestStruct

}  // namespace

TEST(tuple_accessor, tuple_accessor) {

    using namespace detray;

    // detail::make_tuple

    auto s_tuple = detail::make_tuple<std::tuple>(1, 2);

    std::cout << detail::tuple_size<decltype(s_tuple)>::value << std::endl;

    auto t_tuple = detail::make_tuple<thrust::tuple>(1, 2);

    std::cout << detail::tuple_size<decltype(t_tuple)>::value << std::endl;

    // detail::get

    // detail::tuple_size
}
