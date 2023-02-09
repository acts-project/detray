/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

// Detray include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/masks/cuboid3D.hpp"
#include "detray/surface_finders/grid/axis.hpp"
#include "detray/surface_finders/grid/grid.hpp"
#include "detray/surface_finders/grid/populator.hpp"
#include "detray/surface_finders/grid/serializer.hpp"
#include "detray/tools/grid_builder.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <algorithm>
#include <limits>

using namespace detray;
using namespace detray::n_axis;

namespace {

// Algebra definitions
using point3 = __plugin::point3<scalar>;

constexpr scalar inf{std::numeric_limits<scalar>::max()};
constexpr scalar tol{1e-7f};

// Either a data owning or non-owning 3D cartesian multi-axis
template <bool ownership = true, typename containers = host_container_types>
using cartesian_3D = coordinate_axes<cuboid3D::axes<>, ownership, containers>;

// non-owning multi-axis: Takes external containers
bool constexpr is_owning = true;
bool constexpr is_n_owning = false;

// Bin edges for all axes
dvector<scalar> bin_edges = {-10.f, 10.f, -20.f, 20.f, 0.f, 100.f};
// Offsets into edges container and #bins for all axes
dvector<dindex_range> edge_ranges = {{0u, 20u}, {2u, 40u}, {4u, 50u}};

// non-owning multi-axis for the non-owning grid
cartesian_3D<is_n_owning, host_container_types> ax_n_own(&edge_ranges,
                                                         &bin_edges);

// Create some bin data for non-owning grid
template <class populator_t, typename entry_t>
struct bin_content_sequence {

    entry_t entry{0};

    auto operator()() {
        entry += entry_t{1};
        return populator_t::init(entry);
    }
};

/// Test bin content element by element
template <typename grid_t, typename content_t>
void test_content(const grid_t& g, const point3& p, const content_t& expected) {
    dindex i = 0u;
    for (const auto& entry : g.search(p)) {
        ASSERT_NEAR(entry, expected[i++], tol) << " at index " << i - 1u;
    }
}

}  // anonymous namespace

/// Unittest: Test single grid construction
TEST(grid, single_grid) {

    // Owning and non-owning, cartesian, 3-dimensional, replacing grids
    using grid_owning_t =
        grid<cartesian_3D<>, scalar, simple_serializer, replacer>;

    using grid_n_owning_t =
        grid<cartesian_3D<is_n_owning>, scalar, simple_serializer, replacer>;

    using grid_device_t = grid<cartesian_3D<is_owning, device_container_types>,
                               scalar, simple_serializer, replacer>;

    // Fill the bin data for every test
    // bin test entries
    grid_owning_t::bin_storage_type bin_data{};
    bin_data.resize(40'000u);
    std::generate_n(
        bin_data.begin(), 40'000u,
        bin_content_sequence<populator<grid_owning_t::populator_impl>,
                             scalar>());

    // Copy data that will be moved into the data owning types
    dvector<scalar> bin_edges_cp(bin_edges);
    dvector<dindex_range> edge_ranges_cp(edge_ranges);
    grid_owning_t::bin_storage_type bin_data_cp(bin_data);

    // Data-owning axes and grid
    cartesian_3D<is_owning, host_container_types> axes_own(
        std::move(edge_ranges_cp), std::move(bin_edges_cp));
    grid_owning_t grid_own(std::move(bin_data_cp), std::move(axes_own));

    // Check a few basics
    EXPECT_EQ(grid_own.Dim, 3u);
    EXPECT_EQ(grid_own.nbins(), 40'000u);
    auto y_axis = grid_own.get_axis<label::e_y>();
    EXPECT_EQ(y_axis.nbins(), 40u);
    auto z_axis =
        grid_own.get_axis<single_axis<closed<label::e_z>, regular<>>>();
    EXPECT_EQ(z_axis.nbins(), 50u);

    // Create non-owning grid
    grid_n_owning_t grid_n_own(&bin_data, ax_n_own);

    // Test for consistency with owning grid
    EXPECT_EQ(grid_n_own.Dim, grid_own.Dim);
    y_axis = grid_n_own.get_axis<label::e_y>();
    EXPECT_EQ(y_axis.nbins(), grid_own.get_axis<label::e_y>().nbins());
    z_axis = grid_n_own.get_axis<label::e_z>();
    EXPECT_EQ(z_axis.nbins(), grid_own.get_axis<label::e_z>().nbins());

    // Construct a grid from a view
    grid_owning_t::view_type grid_view = get_data(grid_own);
    grid_device_t device_grid(grid_view);

    // Test for consistency with non-owning grid
    EXPECT_EQ(device_grid.Dim, grid_n_own.Dim);
    auto y_axis_dev = device_grid.get_axis<label::e_y>();
    EXPECT_EQ(y_axis_dev.nbins(), grid_n_own.get_axis<label::e_y>().nbins());
    auto z_axis_dev = device_grid.get_axis<label::e_z>();
    EXPECT_EQ(z_axis_dev.nbins(), grid_n_own.get_axis<label::e_z>().nbins());

    // Test const grid view
    /*auto const_grid_view = get_data(const_cast<const
    grid_owning_t&>(grid_own));

    static_assert(
        std::is_same_v<decltype(const_grid_view),
                       typename grid<cartesian_3D<>, const scalar,
                                     simple_serializer, replacer>::view_type>,
        "Const grid view was not correctly constructed!");

    grid<cartesian_3D<is_owning, device_container_types>, const scalar,
         simple_serializer, replacer>
        const_device_grid(const_grid_view);

    static_assert(
        std::is_same_v<typename decltype(const_device_grid)::bin_type,
                       typename replacer::template bin_type<const scalar>>,
        "Const grid was not correctly constructed from view!");*/
}

/// Test bin entry retrieval
TEST(grid, search) {

    // Non-owning, 3D cartesian, replacing grid
    /*using grid_t = grid<scalar, simple_serializer, replacer<>,
                        cartesian_3D, is_n_owning>;
    // init
    grid_t::bin_storage_type bin_data{};
    bin_data.resize(40'000, populator<grid_t::populator_impl>::init<scalar>());

    grid_t grid_3D(&bin_data, ax_n_own);

    axis::multi_bin<3> mbin{2, 1, 0};
    EXPECT_NEAR(*grid_3D.at(mbin), scalar{23}, tol);
    mbin = {14, 3, 23};
    EXPECT_NEAR(*grid_3D.at(mbin), scalar{42}, tol);
    mbin = {0, 16, 7};
    EXPECT_NEAR(*grid_3D.at(mbin), scalar{42}, tol);
    mbin = {6, 31, 44};
    EXPECT_NEAR(*grid_3D.at(mbin), scalar{42}, tol);*/
}

/// Integration test: Test replace population
TEST(grid, replace_population) {

    // Non-owning, 3D cartesian, replacing grid
    using grid_t =
        grid<cartesian_3D<is_n_owning>, scalar, simple_serializer, replacer>;
    // init
    grid_t::bin_storage_type bin_data{};
    bin_data.resize(40'000u, populator<grid_t::populator_impl>::init<scalar>());

    // Create non-owning grid
    grid_t g3r(&bin_data, ax_n_own);

    // Test the initialization
    point3 p = {-10.f, -20.f, 0.f};
    for (int ib0 = 0; ib0 < 20; ++ib0) {
        for (int ib1 = 0; ib1 < 40; ++ib1) {
            for (int ib2 = 0; ib2 < 100; ib2 += 2) {
                p = {static_cast<scalar>(-10 + ib0),
                     static_cast<scalar>(-20 + ib1),
                     static_cast<scalar>(0 + ib2)};
                EXPECT_NEAR(g3r.search(p)[0],
                            std::numeric_limits<scalar>::max(), tol);
            }
        }
    }

    p = {-4.5f, -4.5f, 4.5f};
    // Fill and read
    g3r.populate(p, 3u);
    EXPECT_NEAR(g3r.search(p)[0], static_cast<scalar>(3u), tol);

    // Fill and read two times, fill first 0-99, then 100-199
    for (unsigned int il = 0u; il < 2u; ++il) {
        scalar counter{static_cast<scalar>(il * 100u)};
        for (int ib0 = 0; ib0 < 20; ++ib0) {
            for (int ib1 = 0; ib1 < 40; ++ib1) {
                for (int ib2 = 0; ib2 < 100; ib2 += 2) {
                    p = {static_cast<scalar>(-10 + ib0),
                         static_cast<scalar>(-20 + ib1),
                         static_cast<scalar>(0 + ib2)};
                    g3r.populate(p, counter);
                    EXPECT_NEAR(g3r.search(p)[0], counter, tol);
                    counter += 1.f;
                }
            }
        }
    }

    // A zone test w/o neighbour hood
    /*p = {-4.5, -4.5, 4.5};
    auto test = g2.search(p);
    dvector<dindex> expect = {100u};
    EXPECT_EQ(test, expect);

    // A zone test with neighbour hood
    p = {0.5, 0.5};

    darray<dindex, 2> zone11 = {1u, 1u};
    darray<dindex, 2> zone22 = {2u, 2u};

    test = g2.zone(p, {zone11, zone22}, true);
    expect = {143u, 144u, 145u, 146u, 147u, 153u, 154u, 155u,
              156u, 157u, 163u, 164u, 165u, 166u, 167u};
    EXPECT_EQ(test, expect);

    using grid2cc = grid2<replace_populator, axis::circular, axis::regular,
                          decltype(serializer)>;

    typename grid2cc::axis_p0_type circular{4, -2., 2., host_mr};
    typename grid2cc::axis_p1_type closed{5, 0., 5., host_mr};

    grid2cc g2cc(std::move(circular), std::move(closed), host_mr);
    unsigned int counter = 0;
    for (unsigned icl = 0; icl < 5; ++icl) {
        for (unsigned ici = 0; ici < 4; ++ici) {
            p = {static_cast<scalar>(-1.5 + ici),
                 static_cast<scalar>(0.5 + icl)};
            g2cc.populate(p, counter++);
        }
    }

    // A zone test for circular testing
    p = {1.5, 2.5};
    test = g2cc.zone(p, {zone11, zone11}, true);
    expect = {4u, 6u, 7u, 8u, 10u, 11u, 12u, 14u, 15u};
    EXPECT_EQ(test, expect);*/
}

/// Test bin entry retrieval
TEST(grid, complete_population) {

    // Non-owning, 3D cartesian, completing grid (4 dims and sort)
    using grid_t = grid<cartesian_3D<is_n_owning>, scalar, simple_serializer,
                        completer<4, true>>;
    using bin_content_t = grid_t::bin_type::content_type;

    // init
    grid_t::bin_storage_type bin_data{};
    bin_data.resize(40'000u, populator<grid_t::populator_impl>::init<scalar>());
    // Create non-owning grid
    grid_t g3c(&bin_data, ax_n_own);

    // Test the initialization
    point3 p = {-10.f, -20.f, 0.f};
    grid_t::bin_type invalid =
        populator<grid_t::populator_impl>::init<scalar>();
    for (int ib0 = 0; ib0 < 20; ++ib0) {
        for (int ib1 = 0; ib1 < 40; ++ib1) {
            for (int ib2 = 0; ib2 < 100; ib2 += 2) {
                p = {static_cast<scalar>(-10 + ib0),
                     static_cast<scalar>(-20 + ib1),
                     static_cast<scalar>(0 + ib2)};
                test_content(g3c, p, invalid.content());
            }
        }
    }

    p = {-4.5f, -4.5f, 4.5f};
    bin_content_t expected{4.f, 4.f, 4.f, 4.f};
    // Fill and read
    g3c.populate(p, 4u);
    test_content(g3c, p, expected);

    /*
    auto zone_test = g2.zone(p);
    dvector<dindex> zone_expected = {4u};
    EXPECT_EQ(zone_test, zone_expected);

    g2.populate(p, 2u);
    expected = {4u, 2u, dindex_invalid};
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    g2.populate(p, 7u);
    expected = {4u, 2u, 7u};
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    // Bin is completed, new content is ignored
    g2.populate(p, 16u);
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    darray<dindex, 2> zone00 = {0u, 0u};
    darray<dindex, 2> zone11 = {1u, 1u};

    // Zone test of a complete bin
    zone_test = g2.zone(p, {zone00, zone00});
    zone_expected = {4u, 2u, 7u};
    EXPECT_EQ(zone_test, zone_expected);

    // Fill some other bins
    p = {0.5, -0.5};
    g2.populate(p, 16u);

    p = {0.5, 0.5};
    g2.populate(p, 17u);
    g2.populate(p, 18u);

    zone_test = g2.zone(p, {zone11, zone11});
    zone_expected = {4u, 2u, 7u, 16u, 17u, 18u};
    EXPECT_EQ(zone_test, zone_expected);*/
}

/// Test bin entry retrieval
TEST(grid, regular_attach_population) {

    // Non-owning, 3D cartesian, completing grid (4 dims and sort)
    using grid_t = grid<cartesian_3D<is_n_owning>, scalar, simple_serializer,
                        regular_attacher<4, true>>;
    using bin_content_t = grid_t::bin_type::content_type;

    // init
    grid_t::bin_storage_type bin_data{};
    bin_data.resize(40'000u, populator<grid_t::populator_impl>::init<scalar>());

    // Create non-owning grid
    grid_t g3ra(&bin_data, ax_n_own);

    // Test the initialization
    point3 p = {-10.f, -20.f, 0.f};
    grid_t::bin_type invalid =
        populator<grid_t::populator_impl>::init<scalar>();
    for (int ib0 = 0; ib0 < 20; ++ib0) {
        for (int ib1 = 0; ib1 < 40; ++ib1) {
            for (int ib2 = 0; ib2 < 100; ib2 += 2) {
                p = {static_cast<scalar>(-10 + ib0),
                     static_cast<scalar>(-20 + ib1),
                     static_cast<scalar>(0 + ib2)};
                test_content(g3ra, p, invalid.content());
            }
        }
    }

    p = {-4.5f, -4.5f, 4.5f};
    bin_content_t expected{5.f, inf, inf, inf};
    // Fill and read
    g3ra.populate(p, 5u);
    test_content(g3ra, p, expected);

    /* auto zone_test = g2.zone(p);
     dvector<dindex> zone_expected = {4u};
     EXPECT_EQ(zone_test, zone_expected);

     p = {-0.5, 0.5};
     g2.populate(p, 9u);

     p = {0.5, -0.5};
     g2.populate(p, 1u);

     p = {0.5, 0.5};
     g2.populate(p, 7u);

     expected = {7u};
     test = g2.bin(p);
     EXPECT_EQ(test, expected);

     darray<dindex, 2> zone11 = {1u, 1u};

     zone_test = g2.zone(p, {zone11, zone11}, true);
     zone_expected = {1u, 4u, 7u, 9u};
     EXPECT_EQ(zone_test, zone_expected);*/
}

/*TEST(grids, irregular_replace_population) {

    // Non-owning, 3D cartesian, replacing grid
    using grid_t =
        grid<scalar, simple_serializer, replacer<>, cartesian_3D, is_n_owning>;
    // init
    grid_t::bin_storage_type bin_data{};
    bin_data.resize(40'000, populator<grid_t::populator_impl>::init<scalar>());

    // Create non-owning grid
    grid_t g3r(&bin_data, ax_n_own);
    typename grid2ir::axis_p0_type xaxis{
        {-3, -2., 1, 0.5, 0.7, 0.71, 4., 1000.}, host_mr};
    typename grid2ir::axis_p1_type yaxis{{0.1, 0.8, 0.9, 10., 12., 15.},
                                         host_mr};

    grid2ir g2(std::move(xaxis), std::move(yaxis), host_mr);

    test::point2<detray::scalar> p = {-0.5, 0.5};
    g2.populate(p, 4u);
    EXPECT_EQ(g2.bin(p), 4u);
}*/
