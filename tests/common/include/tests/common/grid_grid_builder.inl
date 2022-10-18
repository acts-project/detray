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
    grid_factory<dindex, simple_serializer, regular_attacher<3>> gf(host_mr);

    auto grid_builder = gf.new_grid_builder();

    // Build from existing mask of a surface
    const scalar minR{0.f};
    const scalar maxR{10.f};
    const scalar minPhi{0.f};
    const scalar maxPhi{M_PI};
    mask<annulus2D<>> ann2{0UL, minR, maxR, minPhi, maxPhi, 0.f, 0.f, 0.f};

    // Grid with correctly initialized axes, but empty bin content
    auto ann_gr = grid_builder.new_grid(ann2, 5, 10);
    const auto& axis_r = ann_gr.template get_axis<label::e_r>();

    // Test axis
    EXPECT_EQ(axis_r.label(), label::e_r);
    EXPECT_EQ(axis_r.shape(), shape::e_open);
    EXPECT_EQ(axis_r.binning(), binning::e_regular);
    EXPECT_EQ(axis_r.nbins(), 5UL);
    EXPECT_NEAR(axis_r.span()[0], 0.f, std::numeric_limits<scalar>::epsilon());
    EXPECT_NEAR(axis_r.span()[1], 10.f, std::numeric_limits<scalar>::epsilon());

    // Test fill a bin to see, if bin content was correctly initialized
    point3 p = {0.5f, 2.f, 0.f};
    vector3 d{};
    auto loc_p = ann_gr.global_to_local(Identity, p, d);
    ann_gr.populate(loc_p, 3UL);

    EXPECT_FLOAT_EQ(ann_gr.search(loc_p)[0], 3UL);


    // Build from parameters
    /*auto cyl_gr = grid_builder.new_cylinder_grid<shape::e_closed, binning::e_regular, binning::e_regular, binning::e_regular>(, 5, 10);
    const auto& axis_r = ann_gr.template get_axis<label::e_r>();

    // Test axis shape
    EXPECT_EQ(axis_r.label(), label::e_r);
    EXPECT_EQ(axis_r.shape(), shape::e_closed);
    EXPECT_EQ(axis_r.binning(), binning::e_regular);

    // N bins
    EXPECT_EQ(axis_r.nbins(), 5UL);
    EXPECT_NEAR(axis_r.span()[0], 0.f, std::numeric_limits<scalar>::epsilon);
    EXPECT_NEAR(axis_r.span()[1], 10.f, std::numeric_limits<scalar>::epsilon);*/
    // Basics
    /*EXPECT_EQ(grid_coll.ngrids(), 3UL);
    EXPECT_EQ(grid_coll.bin_storage().size(), 197UL);
    EXPECT_EQ(grid_coll.axes_storage().size(), 9UL);
    EXPECT_EQ(grid_coll.bin_edges_storage().size(), 18UL);

    // Get a grid instance
    auto single_grid = grid_coll[1];

    static_assert(std::is_same_v<decltype(single_grid), grid_t>,
                  "Grid from collection has wrong type");

    EXPECT_EQ(single_grid.Dim, 3);
    auto r_axis = single_grid.get_axis<label::e_r>();
    EXPECT_EQ(r_axis.nbins(), 1u);
    using z_axis_t = single_axis<open<label::e_z>, regular<>>;
    auto z_axis = single_grid.get_axis<z_axis_t>();
    EXPECT_EQ(z_axis.nbins(), 8u);

    // The generator starts countaing at one instead of zero
    EXPECT_EQ(single_grid.at(0u, 0u, 0u)[0u], 49UL);
    EXPECT_EQ(single_grid.at(0u, 0u, 0u)[1u], inf);
    EXPECT_EQ(single_grid.at(0u, 0u, 0u)[2u], inf);

    auto bin_view = grid_coll[2].at(101u);
    grid_coll[2].populate(101UL, 42UL);
    EXPECT_EQ(bin_view[0u], 102UL + 72UL);
    EXPECT_EQ(bin_view[1u], 42UL);
    EXPECT_EQ(bin_view[2u], inf);*/
}