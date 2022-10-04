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
#include "detray/definitions/containers.hpp"
#include "detray/utils/ranges.hpp"

// Vecmem include(s)
#include "vecmem/containers/device_vector.hpp"
#include "vecmem/containers/jagged_device_vector.hpp"
#include "vecmem/containers/jagged_vector.hpp"

// System include(s)
#include <forward_list>
#include <list>
#include <type_traits>

using namespace detray;

// Test basic ranges definitions
TEST(utils, ranges) {
    //
    // detray containers
    //
    // std::vector
    static_assert(detray::ranges::range_v<dvector<float>>);
    static_assert(not detray::ranges::view<dvector<float>>);
    static_assert(detray::ranges::random_access_range_v<dvector<float>>);
    static_assert(detray::ranges::common_range<dvector<float>>);

    // std::array
    static_assert(detray::ranges::range_v<darray<float, 1>>);
    static_assert(not detray::ranges::view<darray<float, 1>>);
    static_assert(detray::ranges::random_access_range_v<darray<float, 1>>);
    static_assert(detray::ranges::common_range<darray<float, 1>>);

    // std::map
    static_assert(detray::ranges::range_v<dmap<int, float>>);
    static_assert(not detray::ranges::view<dmap<int, float>>);
    static_assert(detray::ranges::bidirectional_range_v<dmap<int, float>>);
    static_assert(not detray::ranges::random_access_range_v<dmap<int, float>>);
    static_assert(detray::ranges::common_range<dmap<int, float>>);

    // std::tuple
    static_assert(not detray::ranges::range_v<dtuple<int, float>>);
    static_assert(not detray::ranges::view<dtuple<int, float>>);

    //
    // vecmem containers
    //

    // vecmem::device_vector
    static_assert(detray::ranges::range_v<vecmem::device_vector<float>>);
    static_assert(not detray::ranges::view<vecmem::device_vector<float>>);
    static_assert(
        detray::ranges::random_access_range_v<vecmem::device_vector<float>>);
    static_assert(detray::ranges::common_range<vecmem::device_vector<float>>);

    // vecmem::jagged_vector
    static_assert(detray::ranges::range_v<djagged_vector<float>>);
    static_assert(not detray::ranges::view<djagged_vector<float>>);
    static_assert(detray::ranges::random_access_range_v<djagged_vector<float>>);
    static_assert(detray::ranges::common_range<djagged_vector<float>>);

    // vecmem::jagged_device_vector
    static_assert(detray::ranges::range_v<vecmem::jagged_device_vector<float>>);
    static_assert(
        not detray::ranges::view<vecmem::jagged_device_vector<float>>);
    static_assert(detray::ranges::random_access_range_v<
                  vecmem::jagged_device_vector<float>>);
    static_assert(
        detray::ranges::common_range<vecmem::jagged_device_vector<float>>);

    //
    // Some additional STL containers
    //

    // std::forward_list
    static_assert(detray::ranges::range_v<std::forward_list<float>>);
    static_assert(not detray::ranges::view<std::forward_list<float>>);
    static_assert(detray::ranges::forward_range_v<std::forward_list<float>>);
    static_assert(
        not detray::ranges::bidirectional_range_v<std::forward_list<float>>);
    static_assert(detray::ranges::common_range<std::forward_list<float>>);

    // std::list
    static_assert(detray::ranges::range_v<std::list<float>>);
    static_assert(not detray::ranges::view<std::list<float>>);
    static_assert(detray::ranges::bidirectional_range_v<std::list<float>>);
    static_assert(not detray::ranges::random_access_range_v<std::list<float>>);
    static_assert(detray::ranges::common_range<std::list<float>>);
}

// Unittest for an empty view
TEST(utils, ranges_empty) {

    auto ev = detray::views::empty<float>();

    // general tests
    static_assert(std::is_copy_assignable_v<decltype(ev)>);
    static_assert(detray::ranges::range_v<decltype(ev)>);
    static_assert(detray::ranges::view<decltype(ev)>);
    static_assert(detray::ranges::random_access_range_v<decltype(ev)>);

    // Test prerequisits for LagacyIterator
    static_assert(
        std::is_copy_constructible_v<typename decltype(ev)::iterator_t>);
    static_assert(std::is_copy_assignable_v<typename decltype(ev)::iterator_t>);
    static_assert(std::is_destructible_v<typename decltype(ev)::iterator_t>);

    // Test inherited member functions
    ASSERT_EQ(ev.size(), 0UL);

    for (const auto i : ev) {
        ASSERT_TRUE(i != i);
    }
}

// Unittest for the generation of a single element sequence
TEST(utils, ranges_single) {

    const dindex value{251UL};

    auto sngl = detray::views::single(value);

    // general tests
    static_assert(std::is_copy_assignable_v<decltype(sngl)>);
    static_assert(detray::ranges::range_v<decltype(sngl)>);
    static_assert(detray::ranges::view<decltype(sngl)>);
    static_assert(detray::ranges::random_access_range_v<decltype(sngl)>);

    // Test prerequisits for LagacyIterator
    static_assert(
        std::is_copy_constructible_v<typename decltype(sngl)::iterator_t>);
    static_assert(
        std::is_copy_assignable_v<typename decltype(sngl)::iterator_t>);
    static_assert(std::is_destructible_v<typename decltype(sngl)::iterator_t>);

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
    static_assert(std::is_copy_assignable_v<decltype(seq)>);
    static_assert(detray::ranges::range_v<decltype(seq)>);
    static_assert(detray::ranges::view<decltype(seq)>);
    static_assert(detray::ranges::input_range_v<decltype(seq)>);
    static_assert(not detray::ranges::forward_range_v<decltype(seq)>);

    // Test prerequisits for LagacyIterator
    static_assert(
        std::is_copy_constructible_v<typename decltype(seq)::iterator_t>);
    static_assert(
        std::is_copy_assignable_v<typename decltype(seq)::iterator_t>);
    static_assert(std::is_destructible_v<typename decltype(seq)::iterator_t>);

    // Test inherited member functions
    ASSERT_EQ(seq.size(), 1UL);

    for (auto& i : seq) {
        check += i;
    }
    ASSERT_EQ(check, single);
}

// Unittest for the generation of a sequence in an interval
TEST(utils, ranges_iota_interval) {

    darray<dindex, 2> interval = {2, 7};

    auto seq = detray::views::iota(interval);

    // general tests
    static_assert(detray::ranges::range_v<decltype(seq)>);
    static_assert(detray::ranges::view<decltype(seq)>);
    static_assert(std::is_copy_assignable_v<decltype(seq)>);
    static_assert(detray::ranges::input_range_v<decltype(seq)>);
    static_assert(not detray::ranges::forward_range_v<decltype(seq)>);

    // Test prerequisits for LagacyIterator
    static_assert(
        std::is_copy_constructible_v<typename decltype(seq)::iterator_t>);
    static_assert(
        std::is_copy_assignable_v<typename decltype(seq)::iterator_t>);
    static_assert(std::is_destructible_v<typename decltype(seq)::iterator_t>);

    // Test inherited member functions
    ASSERT_EQ(seq.size(), 5UL);

    std::vector<dindex> reference = {2, 3, 4, 5, 6};
    std::vector<dindex> check = {};
    for (auto& i : seq) {
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
    static_assert(detray::ranges::range_v<decltype(enumerator)>);
    static_assert(detray::ranges::view<decltype(enumerator)>);
    static_assert(std::is_copy_assignable_v<decltype(enumerator)>);
    static_assert(detray::ranges::random_access_range_v<decltype(enumerator)>);

    // Test prerequisits for LagacyIterator
    static_assert(std::is_copy_constructible_v<
                  typename decltype(enumerator)::iterator_t>);
    static_assert(
        std::is_copy_assignable_v<typename decltype(enumerator)::iterator_t>);
    static_assert(
        std::is_destructible_v<typename decltype(enumerator)::iterator_t>);

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
    static_assert(detray::ranges::range_v<decltype(selected)>);
    static_assert(detray::ranges::view<decltype(selected)>);
    static_assert(std::is_copy_assignable_v<decltype(selected)>);
    static_assert(detray::ranges::random_access_range_v<decltype(selected)>);

    // Test prerequisits for LagacyIterator
    static_assert(
        std::is_copy_constructible_v<typename decltype(selected)::iterator_t>);
    static_assert(
        std::is_copy_assignable_v<typename decltype(selected)::iterator_t>);
    static_assert(
        std::is_destructible_v<typename decltype(selected)::iterator_t>);

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

// Unittest for the joining of multiple ranges
TEST(utils, ranges_join) {

    dvector<dindex> interval_1 = {2, 3, 4};
    dvector<dindex> interval_2 = {7, 8, 9};

    std::vector<dindex> reference = {2, 3, 4, 7, 8, 9};
    std::vector<dindex> check = {};

    auto joined = detray::views::join(interval_1, interval_2);

    // general tests
    static_assert(detray::ranges::range_v<decltype(joined)>);
    static_assert(detray::ranges::view<decltype(joined)>);
    static_assert(std::is_copy_assignable_v<decltype(joined)>);
    static_assert(detray::ranges::random_access_range_v<decltype(joined)>);

    // Test prerequisits for LagacyIterator
    static_assert(
        std::is_copy_constructible_v<typename decltype(joined)::iterator_t>);
    static_assert(
        std::is_copy_assignable_v<typename decltype(joined)::iterator_t>);
    static_assert(
        std::is_destructible_v<typename decltype(joined)::iterator_t>);

    // Test inherited member functions
    ASSERT_EQ(joined[1], 3UL);
    ASSERT_EQ(joined[4], 8UL);
    ASSERT_EQ(joined.size(), 6UL);
    ASSERT_EQ(joined.front(), 2UL);
    ASSERT_EQ(joined.back(), 9UL);

    for (const auto j : joined) {
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

    auto sr = detray::ranges::subrange(seq, interval);

    // general tests
    static_assert(detray::ranges::range_v<decltype(sr)>);
    static_assert(detray::ranges::view<decltype(sr)>);
    static_assert(std::is_copy_assignable_v<decltype(sr)>);
    static_assert(detray::ranges::random_access_range_v<decltype(sr)>);

    // Test prerequisits for LagacyIterator
    static_assert(
        std::is_copy_constructible_v<typename decltype(sr)::iterator_t>);
    static_assert(std::is_copy_assignable_v<typename decltype(sr)::iterator_t>);
    static_assert(std::is_destructible_v<typename decltype(sr)::iterator_t>);

    ASSERT_EQ(sr[1], seq[begin + 1]);
    ASSERT_EQ(sr.size(), 3UL);
    ASSERT_EQ(sr.front(), 1UL);
    ASSERT_EQ(sr.back(), 3UL);

    // non-const iteration
    std::size_t i = 1;
    for (const auto& v : sr) {
        ASSERT_NE(v, 0);
        ASSERT_NE(v, 4);
        ASSERT_EQ(v, seq[i++]);
    }

    // const iteration
    const dvector<int> seq_c(seq);
    i = 1;
    for (const auto& v : detray::ranges::subrange(seq_c, interval)) {
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
    static_assert(detray::ranges::input_range_v<decltype(iota_sr)>);
    static_assert(not detray::ranges::forward_range_v<decltype(iota_sr)>);

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
    static_assert(detray::ranges::random_access_range_v<decltype(enum_sr)>);

    for (const auto [i, v] : enum_sr) {
        ASSERT_EQ(i, v.ui - 1);
    }
}

// Integration test for the picking of indexed elements from another range
TEST(utils, ranges_pick_joined_sequence) {

    darray<dindex, 2> interval_1 = {2, 4};
    darray<dindex, 2> interval_2 = {7, 9};

    // The indices of the iota elements to be picked
    std::vector<dindex> reference = {2, 3, 7, 8};
    std::vector<dindex> check = {};

    struct uint_holder {
        unsigned int ui = 0;
    };

    dvector<uint_holder> seq = {{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}};

    auto indices = detray::views::join(detray::views::iota(interval_1),
                                       detray::views::iota(interval_2));
    auto selected = detray::views::pick(seq, indices);

    // Check iterator category
    static_assert(detray::ranges::input_range_v<decltype(indices)>);
    static_assert(not detray::ranges::forward_range_v<decltype(indices)>);
    static_assert(detray::ranges::input_range_v<decltype(selected)>);
    static_assert(not detray::ranges::forward_range_v<decltype(selected)>);

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