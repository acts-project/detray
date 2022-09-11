/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

// detray test
#include "tests/common/test_defs.hpp"

// detray core
#include "detray/utils/ranges/subrange.hpp"
#include "detray/utils/ranges/views.hpp"

using namespace detray;

// This tests the convenience range_enumeration function: single
/*TEST(utils, views_iota_single) {

    dindex check = 0;
    dindex single = 7;
    for (auto i : detray::ranges::iota_view(single)) {
        check += i;
    }
    ASSERT_EQ(check, single);
}

// This tests the convenience range_enumeration function: range
TEST(utils, views_iota) {

    darray<dindex, 2> range = {2, 7};
    std::vector<dindex> reference = {2, 3, 4, 5, 6, 7};
    std::vector<dindex> check = {};
    for (auto i : detray::ranges::iota_view(range)) {
        check.push_back(i);
    }
    ASSERT_EQ(check, reference);
}*/

// This tests the convenience enumeration function
/*TEST(utils, enumerate) {

    struct uint_holder {
        unsigned int ui = 0;
    };

    dvector<uint_holder> seq = {{0}, {1}, {2}, {3}, {4}, {5}};

    for (auto [i, v] : enumerate(seq)) {
        ASSERT_EQ(i, v.ui);
    }
}*/

// Test the subrange implementation
TEST(utils, ranges_subrange) {

    std::size_t begin = 1;
    std::size_t end = 4;
    std::array<std::size_t, 2> interval{begin, end};

    dvector<int> seq = {0, 1, 2, 3, 4, 5};

    // general tests
    auto sr = detray::ranges::subrange_view(seq, interval);

    ASSERT_TRUE(detray::ranges::range<decltype(sr_c)>::value);
    ASSERT_EQ(sr[1], seq[begin + 1]);
    ASSERT_EQ(sr.size(), 3UL);
    ASSERT_EQ(sr.front(), 1UL);
    ASSERT_EQ(sr.back(), 3UL);

    // non-const iteration
    std::size_t i = 1;
    for (const auto &v : detray::ranges::subrange_view(seq, interval)) {
        ASSERT_NE(v, 4);
        ASSERT_EQ(v, seq[i++]);
    }

    // const iteration
    const dvector<int> seq_c(seq);
    i = 1;
    for (const auto &v : detray::ranges::subrange_view(seq_c, interval)) {
        ASSERT_EQ(v, seq[i++]);
    }
}
