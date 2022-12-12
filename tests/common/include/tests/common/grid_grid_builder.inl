/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

// Detray include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/detectors/detector_metadata.hpp"
#include "detray/masks/masks.hpp"
#include "detray/surface_finders/grid/populator.hpp"
#include "detray/surface_finders/grid/serializer.hpp"
#include "detray/tools/grid_builder.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <limits>

using namespace detray;
using namespace detray::n_axis;

namespace {

using point3 = __plugin::point3<scalar>;
using vector3 = __plugin::vector3<scalar>;

__plugin::transform3<scalar> Identity{};

using test_detector_t = detector<detector_registry::toy_detector>;

}  // anonymous namespace

/// Unittest: Test the construction of a collection of grids
TEST(grid, grid_factory) {

    // Data-owning grid collection
    vecmem::host_memory_resource host_mr;
    auto gr_factory =
        grid_factory<dindex, simple_serializer, regular_attacher<3>>{host_mr};

    // Build from existing mask of a surface
    const scalar minR{0.f};
    const scalar maxR{10.f};
    const scalar minPhi{0.f};
    const scalar maxPhi{M_PI};
    mask<annulus2D<>> ann2{0UL, minR, maxR, minPhi, maxPhi, 0.f, 0.f, 0.f};

    // Grid with correctly initialized axes, but empty bin content
    auto ann_gr = gr_factory.new_grid(ann2, {5, 10});

    // Test axis
    const auto& ann_axis_r = ann_gr.template get_axis<label::e_r>();
    EXPECT_EQ(ann_axis_r.label(), label::e_r);
    EXPECT_EQ(ann_axis_r.bounds(), bounds::e_closed);
    EXPECT_EQ(ann_axis_r.binning(), binning::e_regular);
    EXPECT_EQ(ann_axis_r.nbins(), 5UL);
    EXPECT_NEAR(ann_axis_r.span()[0], 0.f,
                std::numeric_limits<scalar>::epsilon());
    EXPECT_NEAR(ann_axis_r.span()[1], 10.f,
                std::numeric_limits<scalar>::epsilon());

    // Test fill a bin to see, if bin content was correctly initialized
    point3 p = {0.5f, 2.f, 0.f};
    vector3 d{};
    auto loc_p = ann_gr.global_to_local(Identity, p, d);
    ann_gr.populate(loc_p, 3UL);

    EXPECT_FLOAT_EQ(ann_gr.search(loc_p)[0], 3UL);

    // gr_builder.to_string(ann_gr);

    // Build from parameters
    const std::vector<scalar> bin_edges_z{-10.f, -8.f, -6.5f, -1.f,
                                          4.f,   5.f,  6.f,   9.f};
    const std::vector<scalar> bin_edges_phi{};

    auto cyl_gr = gr_factory.template new_grid<cylinder2D<>>(
        {bin_edges_z.front(), bin_edges_z.back(), 0.f,
         2.f * static_cast<scalar>(M_PI)},
        {bin_edges_z.size() - 1, 10UL}, {bin_edges_phi, bin_edges_z},
        std::tuple<circular<label::e_rphi>, closed<label::e_cyl_z>>{},
        std::tuple<regular<>, irregular<>>{});

    // Test axis
    const auto& cyl_axis_z = cyl_gr.template get_axis<label::e_cyl_z>();
    EXPECT_EQ(cyl_axis_z.label(), label::e_cyl_z);
    EXPECT_EQ(cyl_axis_z.bounds(), bounds::e_closed);
    EXPECT_EQ(cyl_axis_z.binning(), binning::e_irregular);
    EXPECT_EQ(cyl_axis_z.nbins(), 7UL);
    EXPECT_NEAR(cyl_axis_z.span()[0], -10.f,
                std::numeric_limits<scalar>::epsilon());
    EXPECT_NEAR(cyl_axis_z.span()[1], 9.f,
                std::numeric_limits<scalar>::epsilon());

    // Test fill a bin to see, if bin content was correctly initialized
    loc_p = cyl_gr.global_to_local(Identity, p, d);
    cyl_gr.populate(loc_p, 33UL);

    EXPECT_FLOAT_EQ(cyl_gr.search(loc_p)[0], 33UL);

    // Build the same cylinder grid from a mask
    const scalar r{5.f};
    const scalar n_half_z{-10.f};
    const scalar p_half_z{9.f};
    mask<cylinder2D<>> cyl2{0UL, r, n_half_z, p_half_z};

    auto cyl_gr2 = gr_factory.template new_grid<circular<label::e_rphi>,
                                                closed<label::e_cyl_z>,
                                                regular<>, irregular<>>(
        cyl2, {bin_edges_z.size() - 1, 10UL}, {bin_edges_phi, bin_edges_z});

    // Get grid collection with the correct allocator
    auto grid_coll = gr_factory.new_collection<decltype(cyl_gr)>();
    EXPECT_TRUE(grid_coll.size() == 0UL);
    grid_coll.push_back(cyl_gr);
    grid_coll.push_back(cyl_gr2);
    EXPECT_TRUE(grid_coll.size() == 2UL);

    EXPECT_FLOAT_EQ(grid_coll[0].search(loc_p)[0], 33UL);
    // gr_factory.to_string(grid_coll[0]);
}

/// Unittest: Test the grid builder
TEST(grid, grid_builder) {

    // cylinder grid type of the toy detector
    using cyl_grid_t =
        grid<coordinate_axes<cylinder2D<>::axes<>, false, host_container_types>,
             test_detector_t::surface_type, simple_serializer,
             regular_attacher<9>>;

    auto gbuilder = grid_builder<test_detector_t, cyl_grid_t,
                                 detray::detail::bin_associator>{};

    // The cylinder portals are at the end of the surface range by construction
    const auto cyl_mask = mask<cylinder2D<>>{0UL, 10.f, -500.f, 500.f};
    std::size_t n_phi_bins{5UL}, n_z_bins{4UL};

    // Build empty grid
    gbuilder.init_grid(cyl_mask, {n_phi_bins, n_z_bins});
    auto cyl_grid = gbuilder();

    const auto& cyl_axis_z = cyl_grid.template get_axis<label::e_cyl_z>();
    EXPECT_EQ(cyl_axis_z.label(), label::e_cyl_z);
    EXPECT_EQ(cyl_axis_z.bounds(), bounds::e_closed);
    EXPECT_EQ(cyl_axis_z.binning(), binning::e_regular);
    EXPECT_EQ(cyl_axis_z.nbins(), 4UL);
    EXPECT_NEAR(cyl_axis_z.span()[0], -500.f,
                std::numeric_limits<scalar>::epsilon());
    EXPECT_NEAR(cyl_axis_z.span()[1], 500.f,
                std::numeric_limits<scalar>::epsilon());
}

/// Integration test: grid builder as volume builder decorator
TEST(grid, decorator_grid_builder) {

    vecmem::host_memory_resource host_mr;

    // cylinder grid type of the toy detector
    using cyl_grid_t =
        grid<coordinate_axes<cylinder2D<>::axes<>, false, host_container_types>,
             test_detector_t::surface_type, simple_serializer,
             regular_attacher<9>>;

    test_detector_t d(host_mr);

    auto vbuilder = std::make_unique<volume_builder<test_detector_t>>();
    auto gbuilder =
        grid_builder<test_detector_t, cyl_grid_t>{std::move(vbuilder)};

    // The cylinder portals are at the end of the surface range by construction
    const auto cyl_mask = mask<cylinder2D<>>{0UL, 10.f, -500.f, 500.f};
    std::size_t n_phi_bins{5UL}, n_z_bins{4UL};

    // Build empty grid
    gbuilder.init_grid(cyl_mask, {n_phi_bins, n_z_bins});

    EXPECT_TRUE(d.volumes().size() == 0);

    // Now init the volume
    gbuilder.init_vol(d, volume_id::e_cylinder,
                      {0., 10., -5., 5., -M_PI, M_PI});

    const auto& vol = d.volumes().back();
    EXPECT_TRUE(d.volumes().size() == 1);
    EXPECT_EQ(vol.index(), 0);
    EXPECT_EQ(vol.id(), volume_id::e_cylinder);

    const auto& cyl_axis_z = gbuilder().template get_axis<label::e_cyl_z>();
    EXPECT_EQ(cyl_axis_z.label(), label::e_cyl_z);
    EXPECT_EQ(cyl_axis_z.nbins(), 4UL);
}