/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

// detray test
#include "tests/common/test_defs.hpp"

// detray core
#include "detray/utils/ranges.hpp"

// System include(s)
#include <type_traits>

using namespace detray;

// Unittest for the generation of a single element sequence
TEST(utils, ranges_single) {

    dindex value{251UL};

    auto sngl = detray::views::single(value);

    // general tests
    ASSERT_TRUE(detray::ranges::range<decltype(sngl)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(sngl)>);
    ASSERT_TRUE(std::is_copy_assignable_v<decltype(sngl)>);

    // Test prerequisits for LagacyIterator
    ASSERT_TRUE(
        std::is_copy_constructible_v<typename decltype(sngl)::iterator_t>);
    ASSERT_TRUE(std::is_copy_assignable_v<typename decltype(sngl)::iterator_t>);
    ASSERT_TRUE(std::is_destructible_v<typename decltype(sngl)::iterator_t>);

    // Check iterator category
    constexpr bool is_correct_itr_tag =
        std::is_same_v<typename decltype(sngl)::iterator_t::iterator_category,
                       std::bidirectional_iterator_tag>;
    ASSERT_TRUE(is_correct_itr_tag);

    // Test inherited member functions
    ASSERT_EQ(sngl[0], value);
    ASSERT_EQ(sngl.size(), 1UL);
    ASSERT_EQ(sngl.front(), 251UL);
    ASSERT_EQ(sngl.back(), 251UL);

    for (auto i : sngl) {
        ASSERT_EQ(251, i);
    }
}

// Unittest for the generation of a single element sequence
TEST(utils, ranges_iota_single) {

    dindex check = 0;
    dindex single = 7;

    auto seq = detray::views::iota(single);

    // general tests
    ASSERT_TRUE(detray::ranges::range<decltype(seq)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(seq)>);
    ASSERT_TRUE(std::is_copy_assignable_v<decltype(seq)>);

    // Test prerequisits for LagacyIterator
    ASSERT_TRUE(
        std::is_copy_constructible_v<typename decltype(seq)::iterator_t>);
    ASSERT_TRUE(std::is_copy_assignable_v<typename decltype(seq)::iterator_t>);
    ASSERT_TRUE(std::is_destructible_v<typename decltype(seq)::iterator_t>);

    // Check iterator category
    constexpr bool is_correct_itr_tag =
        std::is_same_v<typename decltype(seq)::iterator_t::iterator_category,
                       std::input_iterator_tag>;
    ASSERT_TRUE(is_correct_itr_tag);

    // Test inherited member functions
    ASSERT_EQ(seq[1], single + 1);
    ASSERT_EQ(seq.size(), 1UL);
    ASSERT_EQ(seq.front(), 7UL);
    ASSERT_EQ(seq.back(), 7UL);

    for (auto i : seq) {
        check += i;
    }
    ASSERT_EQ(check, single);
}

// Unittest for the generation of a sequence in an interval
TEST(utils, ranges_iota_interval) {

    darray<dindex, 2> interval = {2, 7};

    auto seq = detray::views::iota(interval);

    // general tests
    ASSERT_TRUE(detray::ranges::range<decltype(seq)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(seq)>);
    ASSERT_TRUE(std::is_copy_assignable_v<decltype(seq)>);

    // Test prerequisits for LagacyIterator
    ASSERT_TRUE(
        std::is_copy_constructible_v<typename decltype(seq)::iterator_t>);
    ASSERT_TRUE(std::is_copy_assignable_v<typename decltype(seq)::iterator_t>);
    ASSERT_TRUE(std::is_destructible_v<typename decltype(seq)::iterator_t>);

    // Check iterator category
    constexpr bool is_correct_itr_tag =
        std::is_same_v<typename decltype(seq)::iterator_t::iterator_category,
                       std::input_iterator_tag>;
    ASSERT_TRUE(is_correct_itr_tag);

    // Test inherited member functions
    ASSERT_EQ(seq[1], 3UL);
    ASSERT_EQ(seq.size(), 5UL);
    ASSERT_EQ(seq.front(), 2UL);
    ASSERT_EQ(seq.back(), 6UL);

    std::vector<dindex> reference = {2, 3, 4, 5, 6};
    std::vector<dindex> check = {};
    for (auto i : seq) {
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

    auto enumerator = detray::views::enumerate(seq);

    // general tests
    ASSERT_TRUE(detray::ranges::range<decltype(enumerator)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(enumerator)>);
    ASSERT_TRUE(std::is_copy_assignable_v<decltype(enumerator)>);

    // Test prerequisits for LagacyIterator
    ASSERT_TRUE(std::is_copy_constructible_v<
                typename decltype(enumerator)::iterator_t>);
    ASSERT_TRUE(
        std::is_copy_assignable_v<typename decltype(enumerator)::iterator_t>);
    ASSERT_TRUE(
        std::is_destructible_v<typename decltype(enumerator)::iterator_t>);

    // Check iterator category
    constexpr bool is_correct_itr_tag = std::is_same_v<
        typename decltype(enumerator)::iterator_t::iterator_category,
        std::random_access_iterator_tag>;
    ASSERT_TRUE(is_correct_itr_tag);

    // Test inherited member functions
    const auto [i, v] = enumerator[2];
    ASSERT_EQ(i, 2UL);
    ASSERT_EQ(v.ui, 2UL);
    ASSERT_EQ(enumerator.size(), 6UL);
    const auto [i_front, v_front] = enumerator.front();
    ASSERT_EQ(i_front, 0u);
    ASSERT_EQ(v_front.ui, 0u);
    const auto [i_back, v_back] = enumerator.back();
    ASSERT_EQ(i_back, 5u);
    ASSERT_EQ(v_back.ui, 5u);

    for (auto [j, w] : enumerator) {
        ASSERT_TRUE(j == w.ui);
        ASSERT_EQ(j, w.ui);
    }
}

// Unittest for the chaining of multiple ranges
TEST(utils, ranges_chain) {

    dvector<dindex> interval_1 = {2, 3, 4};
    dvector<dindex> interval_2 = {7, 8, 9};

    std::vector<dindex> reference = {2, 3, 4, 7, 8, 9};
    std::vector<dindex> check = {};

    auto chained = detray::views::chain(interval_1, interval_2);

    // general tests
    ASSERT_TRUE(detray::ranges::range<decltype(chained)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(chained)>);
    ASSERT_TRUE(std::is_copy_assignable_v<decltype(chained)>);

    // Test prerequisits for LagacyIterator
    ASSERT_TRUE(
        std::is_copy_constructible_v<typename decltype(chained)::iterator_t>);
    ASSERT_TRUE(
        std::is_copy_assignable_v<typename decltype(chained)::iterator_t>);
    ASSERT_TRUE(std::is_destructible_v<typename decltype(chained)::iterator_t>);

    // Check iterator category
    constexpr bool is_correct_itr_tag = std::is_same_v<
        typename decltype(chained)::iterator_t::iterator_category,
        std::random_access_iterator_tag>;
    ASSERT_TRUE(is_correct_itr_tag);

    // Test inherited member functions
    ASSERT_EQ(chained[1], 3UL);
    ASSERT_EQ(chained[4], 8UL);
    ASSERT_EQ(chained.size(), 6UL);
    ASSERT_EQ(chained.front(), 2UL);
    ASSERT_EQ(chained.back(), 9UL);

    for (const auto j : chained) {
        check.push_back(j);
    }
    ASSERT_EQ(check.size(), reference.size());
    ASSERT_EQ(check, reference);
}

// Integration test for the picking of indexed elements from another range
TEST(utils, ranges_pick) {

    // The indices of the iota elements to be picked
    std::vector<dindex> indices = {2, 3, 7, 8};
    std::vector<dindex> check = {};

    struct uint_holder {
        unsigned int ui = 0;
    };

    dvector<uint_holder> seq = {{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}};

    auto selected = detray::views::pick(seq, indices);

    // general tests
    ASSERT_TRUE(detray::ranges::range<decltype(selected)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(selected)>);
    ASSERT_TRUE(std::is_copy_assignable_v<decltype(selected)>);

    // Test prerequisits for LagacyIterator
    ASSERT_TRUE(
        std::is_copy_constructible_v<typename decltype(selected)::iterator_t>);
    ASSERT_TRUE(
        std::is_copy_assignable_v<typename decltype(selected)::iterator_t>);
    ASSERT_TRUE(
        std::is_destructible_v<typename decltype(selected)::iterator_t>);

    // Check iterator category
    constexpr bool is_correct_itr_tag = std::is_same_v<
        typename decltype(selected)::iterator_t::iterator_category,
        std::random_access_iterator_tag>;
    ASSERT_TRUE(is_correct_itr_tag);

    // Test inherited member functions
    const auto [i, v] = selected[2];
    ASSERT_EQ(i, 7UL);
    ASSERT_EQ(v.ui, 7UL);
    ASSERT_EQ(selected.size(), 4UL);
    const auto [i_front, v_front] = selected.front();
    ASSERT_EQ(i_front, 2UL);
    ASSERT_EQ(v_front.ui, 2UL);
    const auto [i_back, v_back] = selected.back();
    ASSERT_EQ(i_back, 8UL);
    ASSERT_EQ(v_back.ui, 8UL);

    for (auto [j, w] : selected) {
        ASSERT_TRUE(j == w.ui);
        check.push_back(w.ui);
    }
    ASSERT_EQ(check.size(), indices.size());
    ASSERT_EQ(check, indices);
}

// Unittest for the subrange implementation
TEST(utils, ranges_subrange) {

    std::size_t begin = 1;
    std::size_t end = 4;
    std::array<std::size_t, 2> interval{begin, end};

    dvector<int> seq = {0, 1, 2, 3, 4, 5};

    auto sr = detray::ranges::subrange(seq, interval);

    // general tests
    ASSERT_TRUE(detray::ranges::range<decltype(sr)>::value);
    ASSERT_TRUE(detray::ranges::is_view<decltype(sr)>);
    ASSERT_TRUE(std::is_copy_assignable_v<decltype(sr)>);

    // Test prerequisits for LagacyIterator
    ASSERT_TRUE(
        std::is_copy_constructible_v<typename decltype(sr)::iterator_t>);
    ASSERT_TRUE(std::is_copy_assignable_v<typename decltype(sr)::iterator_t>);
    ASSERT_TRUE(std::is_destructible_v<typename decltype(sr)::iterator_t>);

    // Check iterator category
    constexpr bool is_correct_itr_tag =
        std::is_same_v<typename decltype(sr)::iterator_t::iterator_category,
                       std::random_access_iterator_tag>;
    ASSERT_TRUE(is_correct_itr_tag);

    ASSERT_EQ(sr[1], seq[begin + 1]);
    ASSERT_EQ(sr.size(), 3UL);
    ASSERT_EQ(sr.front(), 1UL);
    ASSERT_EQ(sr.back(), 3UL);

    // non-const iteration
    std::size_t i = 1;
    for (const auto &v : sr) {
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

//
// Integration tests
//

// Integration test for enumeration of a subrange
TEST(utils, ranges_subrange_iota) {

    std::array<std::size_t, 2> seq{1, 10};
    std::array<std::size_t, 2> interval{3, 7};

    auto iota_sr = detray::ranges::subrange(detray::views::iota(seq), interval);

    // Check iterator category
    constexpr bool is_correct_itr_tag = std::is_same_v<
        typename decltype(iota_sr)::iterator_t::iterator_category,
        std::input_iterator_tag>;
    ASSERT_TRUE(is_correct_itr_tag);

    std::size_t i{4};
    for (const auto v : iota_sr) {
        ASSERT_EQ(i++, v);
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

    auto enum_sr =
        detray::views::enumerate(detray::ranges::subrange(seq, interval));

    // Check iterator category
    constexpr bool is_correct_itr_tag = std::is_same_v<
        typename decltype(enum_sr)::iterator_t::iterator_category,
        std::random_access_iterator_tag>;
    ASSERT_TRUE(is_correct_itr_tag);

    for (const auto [i, v] : enum_sr) {
        ASSERT_EQ(i, v.ui - 1);
    }
}

// Integration test for the picking of indexed elements from another range
TEST(utils, ranges_pick_chained_sequence) {

    darray<dindex, 2> interval_1 = {2, 4};
    darray<dindex, 2> interval_2 = {7, 9};

    // The indices of the iota elements to be picked
    std::vector<dindex> reference = {2, 3, 7, 8};
    std::vector<dindex> check = {};

    struct uint_holder {
        unsigned int ui = 0;
    };

    dvector<uint_holder> seq = {{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}};

    auto indices = detray::views::chain(detray::views::iota(interval_1),
                                        detray::views::iota(interval_2));
    auto selected = detray::views::pick(seq, indices);

    // Check iterator category
    bool is_correct_itr_tag = std::is_same_v<
        typename decltype(indices)::iterator_t::iterator_category,
        std::input_iterator_tag>;
    ASSERT_TRUE(is_correct_itr_tag);
    is_correct_itr_tag = std::is_same_v<
        typename decltype(selected)::iterator_t::iterator_category,
        std::input_iterator_tag>;
    ASSERT_TRUE(is_correct_itr_tag);

    // Test inherited member functions
    const auto [i, v] = selected[2];
    ASSERT_EQ(i, 7UL);
    ASSERT_EQ(v.ui, 7UL);
    ASSERT_EQ(selected.size(), 4UL);
    const auto [i_front, v_front] = selected.front();
    ASSERT_EQ(i_front, 2UL);
    ASSERT_EQ(v_front.ui, 2UL);
    const auto [i_back, v_back] = selected.back();
    ASSERT_EQ(i_back, 8UL);
    ASSERT_EQ(v_back.ui, 8UL);

    for (auto [j, w] : selected) {
        ASSERT_TRUE(j == w.ui);
        check.push_back(w.ui);
    }
    ASSERT_EQ(check.size(), reference.size());
    ASSERT_EQ(check, reference);
}