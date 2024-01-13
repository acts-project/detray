/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/masks/cylinder3D.hpp"
#include "detray/surface_finders/grid/axis.hpp"
#include "detray/surface_finders/grid/grid.hpp"
#include "detray/surface_finders/grid/grid_collection.hpp"
#include "detray/surface_finders/grid/populators.hpp"
#include "detray/surface_finders/grid/serializers.hpp"
#include "detray/test/types.hpp"
#include "detray/tools/grid_builder.hpp"

// System include(s)
#include <algorithm>
#include <climits>

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;
using namespace detray::n_axis;

namespace {

// non-owning multi-axis: Takes external containers
bool constexpr is_n_owning = false;

constexpr dindex inf{std::numeric_limits<dindex>::max()};

// Create some bin data for non-owning grid
template <typename populator_t, typename bin_t>
struct bin_content_sequence {
    using entry_t = typename bin_t::entry_type;
    entry_t entry{0};

    auto operator()() {
        entry += entry_t{1};
        bin_t bin{};
        populator_t{}(bin, entry);
        return bin;
    }
};

}  // anonymous namespace

/// Unittest: Test the construction of a collection of grids
GTEST_TEST(detray_grid, grid_collection) {

    // grid type

    // Non-owning grid type with array<dindex, 3> as bin content
    using cylindrical_3D = coordinate_axes<cylinder3D::axes<>, is_n_owning>;
    using grid_t = grid<cylindrical_3D, bins::static_array<dindex, 3>>;

    // Build test data

    // Offsets into edges container and #bins for all axes
    dvector<dindex_range> edge_ranges = {{0u, 2u},  {2u, 4u},  {4u, 6u},
                                         {6u, 1u},  {8u, 3u},  {10u, 8u},
                                         {12u, 5u}, {14u, 5u}, {16u, 5u}};

    // Bin edges for all axes
    dvector<scalar> bin_edges = {-10, 10., -20., 20., 0.,  120., -5., 5., -15.,
                                 15., 0.,  50.,  -15, 15., -35., 35., 0., 550.};

    // Bin test entries
    grid_t::bin_container_type bin_data{};
    bin_data.resize(197u);
    std::generate_n(
        bin_data.begin(), 197u,
        bin_content_sequence<attach<>, typename grid_t::bin_type>());
    dvector<dindex> grid_offsets = {0u, 48u, 72u};

    // Data-owning grid collection
    auto grid_coll =
        grid_collection<grid_t>(std::move(grid_offsets), std::move(bin_data),
                                std::move(edge_ranges), std::move(bin_edges));

    // Tests

    // Basics
    EXPECT_EQ(grid_coll.size(), 3u);
    EXPECT_EQ(grid_coll.bin_storage().size(), 197u);
    EXPECT_EQ(grid_coll.axes_storage().size(), 9u);
    EXPECT_EQ(grid_coll.bin_edges_storage().size(), 18u);

    // Get a grid instance
    auto single_grid = grid_coll[1];

    static_assert(std::is_same_v<decltype(single_grid), grid_t>,
                  "Grid from collection has wrong type");

    EXPECT_EQ(single_grid.Dim, 3);
    EXPECT_EQ(single_grid.nbins(), 24u);
    auto r_axis = single_grid.get_axis<label::e_r>();
    EXPECT_EQ(r_axis.nbins(), 1u);
    auto phi_axis = single_grid.get_axis<label::e_phi>();
    EXPECT_EQ(phi_axis.nbins(), 3u);
    using z_axis_t = single_axis<closed<label::e_z>, regular<>>;
    auto z_axis = single_grid.get_axis<z_axis_t>();
    EXPECT_EQ(z_axis.nbins(), 8u);

    // The generator starts counting at one instead of zero
    EXPECT_EQ(single_grid.bin(0u, 0u, 0u)[0u], 49u);
    EXPECT_EQ(single_grid.bin(0u, 0u, 0u)[1u], inf);
    EXPECT_EQ(single_grid.bin(0u, 0u, 0u)[2u], inf);

    // Test the bin view
    auto& bin_view = grid_coll[2].bin(101u);
    grid_coll[2].template populate<attach<>>(101u, 42u);
    EXPECT_EQ(bin_view[0u], 102u + 72u);
    EXPECT_EQ(bin_view[1u], 42u);
    EXPECT_EQ(bin_view[2u], inf);

    // Test the global bin iteration. Take the middle grid!
    auto seq = detray::views::iota(49, 73);
    auto flat_bin_view = single_grid.all();
    EXPECT_EQ(seq.size(), 24u);
    EXPECT_EQ(flat_bin_view.size(), 24u);
    EXPECT_TRUE(
        std::equal(flat_bin_view.begin(), flat_bin_view.end(), seq.begin()));

    auto grid_coll_view = get_data(grid_coll);
    static_assert(std::is_same_v<decltype(grid_coll_view),
                                 typename grid_collection<grid_t>::view_type>,
                  "Grid collection view incorrectly assembled");

    const grid_collection<grid_t>& const_coll = grid_coll;
    auto const_coll_view = get_data(const_coll);
    static_assert(
        std::is_same_v<decltype(const_coll_view),
                       typename grid_collection<grid_t>::const_view_type>,
        "Grid collection const view incorrectly assembled");
}
