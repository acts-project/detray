/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <climits>

// detray core
#include "detray/definitions/indexing.hpp"
#include "detray/surface_finders/grid/populator.hpp"
#include "detray/test/types.hpp"

using namespace detray;

namespace {

constexpr dindex inf{std::numeric_limits<dindex>::max()};

template <class populator_t, typename entry_t>
struct increment {
    using bin_t = typename populator_t::template bin_type<entry_t>;

    bin_t stored_content;
    entry_t entry;

    increment() : entry{0} {}
    bin_t operator()() {
        entry += entry_t{1};
        return populator_t::init(entry);
    }
};

/// Test bin content entry by entry
template <typename populator_t, typename storage_t, typename content_t>
void test_content(populator_t& p, storage_t& storage, dindex bin_idx,
                  content_t& expected) {
    dindex i = 0;
    for (const auto& entry : p.view(storage, bin_idx)) {
        EXPECT_EQ(entry, expected[i++]);
    }
}

}  // anonymous namespace

/// Replace populator
GTEST_TEST(detray_grid, replace_populator) {
    // No sorting, dindex entries in backend storage, vecmem::vector
    using populator_t = populator<replacer>;
    populator_t replace_populator;

    // Create some bin data
    dvector<populator_t::template bin_type<dindex>> bin_data{};
    bin_data.reserve(50);
    std::generate_n(bin_data.begin(), 50, increment<populator_t, dindex>());

    // Check test setup
    EXPECT_EQ(bin_data[2].content(), 3u);
    EXPECT_EQ(replace_populator.view(bin_data, 2)[0], 3u);
    EXPECT_EQ(bin_data[42].content(), 43u);
    EXPECT_EQ(replace_populator.view(bin_data, 42)[0], 43u);

    // Replace some bin entries
    dindex entry{15};

    replace_populator(bin_data, 2, entry);
    EXPECT_EQ(bin_data[2].content(), 15u);
    EXPECT_EQ(replace_populator.view(bin_data, 2)[0], 15u);

    replace_populator(bin_data, 2, entry);
    EXPECT_EQ(bin_data[2].content(), 15u);
    EXPECT_EQ(replace_populator.view(bin_data, 2)[0], 15u);
}

/// Complete populator
GTEST_TEST(detray_grid, complete_populator) {

    // No sorting, 4 dims, dindex entries in backend storage, std::array
    using populator_t = populator<completer<4>>;
    populator_t complete_populator{};

    // Create some bin data
    dvector<populator_t::template bin_type<dindex>> bin_data{};
    bin_data.reserve(50);
    std::generate_n(bin_data.begin(), 50, increment<populator_t, dindex>());

    // Check test setup
    populator_t::template bin_type<dindex>::content_type stored = {3u, inf, inf,
                                                                   inf};
    EXPECT_EQ(bin_data[2].content(), stored);
    test_content(complete_populator, bin_data, 2, stored);
    stored = {43u, inf, inf, inf};
    EXPECT_EQ(bin_data[42].content(), stored);
    test_content(complete_populator, bin_data, 42, stored);

    // Fill up some bin entries
    dindex entry{15};

    stored = {3u, 15u, 15u, 15u};
    complete_populator(bin_data, 2, entry);
    EXPECT_EQ(bin_data[2].content(), stored);
    test_content(complete_populator, bin_data, 2, stored);

    stored = {43u, 15u, 15u, 15u};
    complete_populator(bin_data, 42, entry);
    EXPECT_EQ(bin_data[42].content(), stored);
    test_content(complete_populator, bin_data, 42, stored);

    // Do sorting, 4 dims, dindex entries in backend storage, std::array
    populator<completer<4, true>> sort_complete_populator;

    stored = {2u, 15u, 15u, 15u};
    sort_complete_populator(bin_data, 1, entry);
    EXPECT_EQ(bin_data[1].content(), stored);
    test_content(sort_complete_populator, bin_data, 1, stored);

    stored = {15u, 15u, 15u, 42u};
    sort_complete_populator(bin_data, 41, entry);
    EXPECT_EQ(bin_data[41].content(), stored);
    test_content(sort_complete_populator, bin_data, 41, stored);
}

/// Regular attach populator
GTEST_TEST(detray_grid, regular_attach_populator) {

    // No sorting, 4 dims, dindex entries in backend storage, std::array
    using populator_t = populator<regular_attacher<4>>;
    populator_t reg_attach_populator{};

    // Create some bin data
    dvector<populator_t::template bin_type<dindex>> bin_data{};
    bin_data.reserve(50);
    std::generate_n(bin_data.begin(), 50, increment<populator_t, dindex>());

    // Check test setup
    populator_t::template bin_type<dindex>::content_type stored = {3u, inf, inf,
                                                                   inf};
    EXPECT_EQ(bin_data[2].content(), stored);
    test_content(reg_attach_populator, bin_data, 2, stored);
    stored = {43u, inf, inf, inf};
    EXPECT_EQ(bin_data[42].content(), stored);
    test_content(reg_attach_populator, bin_data, 42, stored);

    // Attach some bin entries
    dindex entry1{15}, entry2{8};

    stored = {3u, 15u, 8u, inf};
    reg_attach_populator(bin_data, 2, entry1);
    reg_attach_populator(bin_data, 2, entry2);
    EXPECT_EQ(bin_data[2].content(), stored);
    test_content(reg_attach_populator, bin_data, 2, stored);

    stored = {43u, 15u, 8u, 16u};
    reg_attach_populator(bin_data, 42, entry1);
    reg_attach_populator(bin_data, 42, entry2);
    reg_attach_populator(bin_data, 42, entry1 + 1);
    EXPECT_EQ(bin_data[42].content(), stored);
    test_content(reg_attach_populator, bin_data, 42, stored);

    // Do sorting, 4 dims, dindex entries in backend storage, std::array
    populator<regular_attacher<4, true>> sort_reg_attach_populator;

    stored = {2u, 8u, 9u, 15u};
    sort_reg_attach_populator(bin_data, 1, entry1);
    sort_reg_attach_populator(bin_data, 1, entry2);
    sort_reg_attach_populator(bin_data, 1, entry2 + 1);
    EXPECT_EQ(bin_data[1].content(), stored);
    test_content(sort_reg_attach_populator, bin_data, 1, stored);

    stored = {8u, 15u, 42u, inf};
    sort_reg_attach_populator(bin_data, 41, entry1);
    sort_reg_attach_populator(bin_data, 41, entry2);
    EXPECT_EQ(bin_data[41].content(), stored);
    test_content(sort_reg_attach_populator, bin_data, 41, stored);
}
