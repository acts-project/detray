/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

// Detray include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/masks/cylinder3D.hpp"
#include "detray/surface_finders/grid/axis.hpp"
#include "detray/surface_finders/grid/grid.hpp"
#include "detray/surface_finders/grid/grid_builder.hpp"
#include "detray/surface_finders/grid/grid_collection.hpp"
#include "detray/surface_finders/grid/populator.hpp"
#include "detray/surface_finders/grid/serializer.hpp"

// System include(s)
#include <algorithm>
#include <climits>

using namespace detray;

namespace {

// non-owning multi-axis: Takes external containers
bool constexpr is_owning = true;
bool constexpr is_n_owning = false;

// Create some bin data for non-owning grid
template <class populator_t, typename entry_t>
struct increment {
    using bin_t = typename populator_t::template bin_type<entry_t>;

    bin_t stored_entry;
    entry_t entry;

    increment() : entry{0} {}
    bin_t operator()() {
        entry += entry_t{1};
        return populator_t::init(entry);
    }
};

/// Test bin content element by element
/*template <typename grid_t, typename content_t>
void test_content(const grid_t& g, const point3& p, const content_t& expected) {
    dindex i = 0;
    for (const auto& entry : g.search(p)) {
        ASSERT_FLOAT_EQ(entry, expected[i++]);
    }
}*/

}  // anonymous namespace

/// Unittest: Test the construction of a collection of grids
TEST(grid, grid_collection) {

    // grid type

    // Non-owning grid type with array<dindex, 3> as bin content
    using cylindrical_3D =
        coordinate_axes<cylinder3D::axes<>, is_n_owning, host_container_types>;
    using grid_t =
        grid<cylindrical_3D, dindex, simple_serializer, regular_attacher<3>>;

    // Build test data

    // Offsets into edges container and #bins for all axes
    dvector<dindex_range> edge_ranges = {{0u, 20u},  {2u, 40u},  {4u, 60u},
                                         {6u, 10u},  {8u, 30u},  {10u, 50u},
                                         {12u, 15u}, {14u, 35u}, {16u, 55u}};

    // Bin edges for all axes
    dvector<scalar> bin_edges = {-10, 10., -20., 20., 0.,  120., -5., 5., -15.,
                                 15., 0.,  50.,  -15, 15., -35., 35., 0., 550.};

    // Bin test entries
    grid_t::bin_storage_type bin_data{};
    bin_data.resize(315);
    std::generate_n(bin_data.begin(), 315u,
                    increment<grid_t::populator_type, dindex>());

    // Data-owning grid collection
    auto gr_coll = grid_collection<grid_t>(
        std::move(bin_data), std::move(edge_ranges), std::move(bin_edges));

    // Tests

    // Basics
    EXPECT_EQ(gr_coll.bin_data().size(), 315UL);
    EXPECT_EQ(gr_coll.axes_data().size(), 9UL);
    EXPECT_EQ(gr_coll.bin_edges_data().size(), 18UL);
}