/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

// Detray include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/masks/masks.hpp"
#include "detray/surface_finders/grid/grid_builder.hpp"
#include "detray/surface_finders/grid/populator.hpp"
#include "detray/surface_finders/grid/serializer.hpp"

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

}  // anonymous namespace

/// Unittest: Test the construction of a collection of grids
TEST(grid, grid_builder) {

    // Data-owning grid collection
    vecmem::host_memory_resource host_mr;
    grid_builder<dindex, simple_serializer, regular_attacher<3>> gr_builder(
        host_mr);

    auto grid_factory = gr_builder.new_factory();

    // Build from existing mask of a surface
    const scalar minR{0.f};
    const scalar maxR{10.f};
    const scalar minPhi{0.f};
    const scalar maxPhi{M_PI};
    mask<annulus2D<>> ann2{0UL, minR, maxR, minPhi, maxPhi, 0.f, 0.f, 0.f};

    // Grid with correctly initialized axes, but empty bin content
    auto ann_gr = grid_factory.new_grid(ann2, 5, 10);

    // Test axis
    const auto& ann_axis_r = ann_gr.template get_axis<label::e_r>();
    EXPECT_EQ(ann_axis_r.label(), label::e_r);
    EXPECT_EQ(ann_axis_r.shape(), shape::e_open);
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

    auto cyl_gr = grid_factory.template new_grid<cylinder2D<>, shape::e_closed,
                                                 regular, irregular>(
        {bin_edges_z.front(), bin_edges_z.back(), 0.f,
         2.f * static_cast<scalar>(M_PI)},
        {bin_edges_z.size() - 1, 10UL}, {bin_edges_phi, bin_edges_z});

    // Test axis
    const auto& cyl_axis_z = cyl_gr.template get_axis<label::e_cyl_z>();
    EXPECT_EQ(cyl_axis_z.label(), label::e_cyl_z);
    EXPECT_EQ(cyl_axis_z.shape(), shape::e_closed);
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

    auto cyl_gr2 =
        grid_factory.template new_grid<shape::e_closed, regular, irregular>(
            cyl2, bin_edges_z.size() - 1, 10UL, bin_edges_phi, bin_edges_z);

    // Get grid collection with the correct allocator
    auto grid_coll = gr_builder.new_collection<decltype(cyl_gr)>();
    EXPECT_TRUE(grid_coll.size() == 0UL);
    grid_coll.push_back(cyl_gr);
    grid_coll.push_back(cyl_gr2);
    EXPECT_TRUE(grid_coll.size() == 2UL);

    EXPECT_FLOAT_EQ(grid_coll[0].search(loc_p)[0], 33UL);
    // gr_builder.to_string(grid_coll[0]);
}