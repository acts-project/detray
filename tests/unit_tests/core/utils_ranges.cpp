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
#include "detray/utils/ranges.hpp"

using namespace detray;

// Unittest for the generation of a single element sequence
TEST(utils, ranges_single) {

    dindex value{251UL};

    // general tests
    auto sngl = detray::views::single(value);
    ASSERT_TRUE(detray::ranges::range<decltype(sngl)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(sngl)>);
    ASSERT_EQ(sngl[0], value);
    ASSERT_EQ(sngl.size(), 1UL);
    ASSERT_EQ(sngl.front(), 251UL);
    ASSERT_EQ(sngl.back(), 251UL);

    for (auto i : detray::views::single(value)) {
        ASSERT_EQ(251, i);
    }
}

// Unittest for the generation of a single element sequence
TEST(utils, ranges_iota_single) {

    dindex check = 0;
    dindex single = 7;

    // general tests
    auto seq = detray::views::iota(single);
    ASSERT_TRUE(detray::ranges::range<decltype(seq)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(seq)>);
    ASSERT_EQ(seq[1], single + 1);
    ASSERT_EQ(seq.size(), 1UL);
    ASSERT_EQ(seq.front(), 7UL);

    for (auto i : detray::views::iota(single)) {
        check += i;
    }
    ASSERT_EQ(check, single);
}

// Unittest for the generation of a sequence in an interval
TEST(utils, ranges_iota_interval) {

    darray<dindex, 2> interval = {2, 7};

    // general tests
    auto seq = detray::views::iota(interval);
    ASSERT_TRUE(detray::ranges::range<decltype(seq)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(seq)>);
    ASSERT_EQ(seq[1], 3UL);
    ASSERT_EQ(seq.size(), 5UL);
    ASSERT_EQ(seq.front(), 2UL);

    std::vector<dindex> reference = {2, 3, 4, 5, 6};
    std::vector<dindex> check = {};
    for (auto i : detray::views::iota(interval)) {
        check.push_back(i);
    }
    ASSERT_EQ(check.size(), reference.size());
    ASSERT_EQ(check, reference);
}

// Unittest for the convenience enumeration of a range
TEST(utils, ranges_enumerate) {

    struct uint_holder {
        unsigned int ui = 0;
    };

    dvector<uint_holder> seq = {{0}, {1}, {2}, {3}, {4}, {5}};

    // general tests
    auto enumerator = detray::views::enumerate(seq);
    ASSERT_TRUE(detray::ranges::range<decltype(enumerator)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(enumerator)>);
    ASSERT_EQ(enumerator.size(), 6UL);
    const auto [i_front, v_front] = enumerator.front();
    ASSERT_EQ(i_front, 0u);
    ASSERT_EQ(v_front.ui, 0u);

    for (auto [i, v] : detray::views::enumerate(seq)) {
        ASSERT_EQ(i, v.ui);
    }
}

// Unittest for the chaining of multiple ranges
TEST(utils, ranges_chain) {

    darray<dindex, 2> interval_1 = {2, 5};
    darray<dindex, 2> interval_2 = {7, 10};

    std::vector<dindex> reference = {2, 3, 4, 7, 8, 9};
    std::vector<dindex> check = {};
    // general tests
    auto chained = detray::views::chain(detray::views::iota(interval_1),
                                        detray::views::iota(interval_2));
    ASSERT_TRUE(detray::ranges::range<decltype(chained)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(chained)>);
    ASSERT_EQ(chained.front(), 2UL);

    for (const auto j : chained) {
        check.push_back(j);
    }
    ASSERT_EQ(check.size(), reference.size());
    ASSERT_EQ(check, reference);
}

// Unittest for the subrange implementation
TEST(utils, ranges_subrange) {

    std::size_t begin = 1;
    std::size_t end = 4;
    std::array<std::size_t, 2> interval{begin, end};

    dvector<int> seq = {0, 1, 2, 3, 4, 5};

    // general tests
    auto sr = detray::ranges::subrange(seq, interval);

    ASSERT_TRUE(detray::ranges::range<decltype(sr)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(sr)>);
    ASSERT_EQ(sr[1], seq[begin + 1]);
    ASSERT_EQ(sr.size(), 3UL);
    ASSERT_EQ(sr.front(), 1UL);
    ASSERT_EQ(sr.back(), 3UL);

    // non-const iteration
    std::size_t i = 1;
    for (const auto &v : detray::ranges::subrange(seq, interval)) {
        ASSERT_NE(v, 0);
        ASSERT_NE(v, 4);
        ASSERT_EQ(v, seq[i++]);
    }

    // const iteration
    const dvector<int> seq_c(seq);
    i = 1;
    for (const auto &v : detray::ranges::subrange(seq_c, interval)) {
        ASSERT_EQ(v, seq[i++]);
    }
}

// Integration test for enumeration of a subrange
TEST(utils, ranges_enumerated_subrange) {

    struct uint_holder {
        unsigned int ui = 0;
    };

    dvector<uint_holder> seq = {{0}, {1}, {2}, {3}, {4}, {5}};

    std::size_t begin = 1;
    std::size_t end = 4;
    std::array<std::size_t, 2> interval{begin, end};

    for (const auto [i, v] :
         detray::views::enumerate(detray::ranges::subrange(seq, interval))) {
        ASSERT_EQ(i, v.ui - 1);
    }
}