/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <climits>

// detray test
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// detray core
#include "detray/definitions/indexing.hpp"
#include "detray/masks/masks.hpp"
#include "detray/surface_finders/grid/axis.hpp"

using namespace detray;
using namespace detray::n_axis;

namespace {

using point3 = __plugin::point3<scalar>;

// Alias for testing
template <bool ownership, typename containers>
using cartesian_3D =
    multi_axis<ownership, cartesian3<__plugin::transform3<scalar>>,
               single_axis<closed<label::e_x>, regular<containers, scalar>>,
               single_axis<closed<label::e_y>, regular<containers, scalar>>,
               single_axis<closed<label::e_z>, regular<containers, scalar>>>;

// Floating point comparison
constexpr scalar tol{1e-5};

}  // anonymous namespace

TEST(grid, open_regular_axis) {

    // Lower bin edges: min and max bin edge for the regular axis
    vecmem::vector<scalar> bin_edges = {-10, -5, -3, 7, 7, 14, 20};
    // Regular axis: first entry is the offset in the bin edges vector (2), the
    // second entry is the number of bins (10): Lower and upper bin edges of
    // the max and min bin are -3 and 7 => stepsize is (7 - (-3)) / 10 = 1
    // ... -3, -2, -1,  0,  1,  2,  3,  4,  5,  6,  7 ...
    //  [0]  [1] [2] [3] [4] [5] [6] [7] [8] [9] [10] [11]
    dindex_range edge_range = {2, 10};

    // An open regular x-axis
    single_axis<open<label::e_x>, regular<>> or_axis{&edge_range, &bin_edges};

    // Test axis shape
    EXPECT_EQ(or_axis.label(), n_axis::label::e_x);
    EXPECT_EQ(or_axis.shape(), n_axis::shape::e_open);
    EXPECT_EQ(or_axis.binning(), n_axis::binning::e_regular);

    // N bins
    EXPECT_EQ(or_axis.nbins(), 10u);
    EXPECT_NEAR(or_axis.span()[0], scalar{-3}, tol);
    EXPECT_NEAR(or_axis.span()[1], scalar{7}, tol);
    // Axis bin access
    // Axis is open: Every value smaller than -3 is mapped into bin 0:
    // which is (inf, -3)
    EXPECT_EQ(or_axis.bin(-4.), 0u);
    EXPECT_EQ(or_axis.bin(2.5), 6u);
    // Axis is open: Every value greater than 7 is mapped into bin 11:
    // which is (7, inf)
    EXPECT_EQ(or_axis.bin(8.), 11u);

    // Axis range access - binned (symmetric & asymmetric)
    darray<dindex, 2> nhood00i = {0u, 0u};
    darray<dindex, 2> nhood01i = {0u, 1u};
    darray<dindex, 2> nhood11i = {1u, 1u};
    darray<dindex, 2> nhood44i = {4u, 4u};
    darray<dindex, 2> nhood55i = {5u, 5u};

    dindex_range expected_range = {6u, 6u};
    EXPECT_EQ(or_axis.range(2.5, nhood00i), expected_range);
    expected_range = {6u, 7u};
    EXPECT_EQ(or_axis.range(2.5, nhood01i), expected_range);
    expected_range = {5u, 7u};
    EXPECT_EQ(or_axis.range(2.5, nhood11i), expected_range);
    expected_range = {2u, 10u};
    EXPECT_EQ(or_axis.range(2.5, nhood44i), expected_range);
    expected_range = {1u, 11u};
    EXPECT_EQ(or_axis.range(2.5, nhood55i), expected_range);
    expected_range = {1u, 9u};
    EXPECT_EQ(or_axis.range(1.5, nhood44i), expected_range);
    expected_range = {4u, 11u};
    EXPECT_EQ(or_axis.range(5.5, nhood55i), expected_range);

    // Axis range access - scalar (symmteric & asymmetric)
    darray<scalar, 2> nhood00s = {0., 0.};
    darray<scalar, 2> epsilon = {0.01, 0.01};
    darray<scalar, 2> nhood11s = {1., 1.};
    darray<scalar, 2> nhoodAlls = {10., 10.};

    expected_range = {6u, 6u};
    EXPECT_EQ(or_axis.range(2.5, nhood00s), expected_range);
    EXPECT_EQ(or_axis.range(2.5, epsilon), expected_range);
    expected_range = {5u, 7u};
    EXPECT_EQ(or_axis.range(2.5, nhood11s), expected_range);
    expected_range = {0u, 11u};
    EXPECT_EQ(or_axis.range(2.5, nhoodAlls), expected_range);
}

TEST(grid, closed_regular_axis) {

    // Lower bin edges: min and max bin edge for the regular axis
    vecmem::vector<scalar> bin_edges = {-10, -3, -3, 7, 7, 14};
    // Regular axis: first entry is the offset in the bin edges vector (2), the
    // second entry is the number of bins (10): Lower and upper bin edges of
    // the max and min bin are -3 and 7 => stepsize is (7 - (-3)) / 10 = 1
    // -3, -2, -1,  0,  1,  2,  3,  4,  5,  6,  7
    //   [0] [1] [2] [3] [4] [5] [6] [7] [8] [9]
    dindex_range edge_range = {2, 10};

    // A closed regular r-axis
    single_axis<closed<label::e_r>, regular<>> cr_axis{&edge_range, &bin_edges};

    // Test axis shape
    EXPECT_EQ(cr_axis.label(), n_axis::label::e_r);
    EXPECT_EQ(cr_axis.shape(), n_axis::shape::e_closed);
    EXPECT_EQ(cr_axis.binning(), n_axis::binning::e_regular);

    // N bins
    EXPECT_EQ(cr_axis.nbins(), 10u);
    EXPECT_NEAR(cr_axis.span()[0], scalar{-3}, tol);
    EXPECT_NEAR(cr_axis.span()[1], scalar{7}, tol);
    // Axis bin access
    // Axis is closed: Every value smaller than -3 is mapped into bin 0:
    // which is [-3, -2)
    EXPECT_EQ(cr_axis.bin(-4.), 0u);
    EXPECT_EQ(cr_axis.bin(2.5), 5u);
    // Axis is closed: Every value greater than 7 is mapped into bin 9:
    // which is (6, 7]
    EXPECT_EQ(cr_axis.bin(8.), 9u);

    // Axis range access - binned (symmetric & asymmetric)
    darray<dindex, 2> nhood00i = {0u, 0u};
    darray<dindex, 2> nhood01i = {0u, 1u};
    darray<dindex, 2> nhood11i = {1u, 1u};
    darray<dindex, 2> nhood44i = {4u, 4u};
    darray<dindex, 2> nhood55i = {5u, 5u};

    dindex_range expected_range = {5u, 5u};
    EXPECT_EQ(cr_axis.range(2.5, nhood00i), expected_range);
    expected_range = {5u, 6u};
    EXPECT_EQ(cr_axis.range(2.5, nhood01i), expected_range);
    expected_range = {4u, 6u};
    EXPECT_EQ(cr_axis.range(2.5, nhood11i), expected_range);
    expected_range = {1u, 9u};
    EXPECT_EQ(cr_axis.range(2.5, nhood44i), expected_range);
    expected_range = {0u, 9u};
    EXPECT_EQ(cr_axis.range(2.5, nhood55i), expected_range);
    expected_range = {0u, 8u};
    EXPECT_EQ(cr_axis.range(1.5, nhood44i), expected_range);
    expected_range = {3u, 9u};
    EXPECT_EQ(cr_axis.range(5.5, nhood55i), expected_range);

    // Axis range access - scalar (symmteric & asymmetric)
    darray<scalar, 2> nhood00s = {0., 0.};
    darray<scalar, 2> epsilon = {0.01, 0.01};
    darray<scalar, 2> nhood11s = {1., 1.};
    darray<scalar, 2> nhoodAlls = {10., 10.};

    expected_range = {5u, 5u};
    EXPECT_EQ(cr_axis.range(2.5, nhood00s), expected_range);
    EXPECT_EQ(cr_axis.range(2.5, epsilon), expected_range);
    expected_range = {4u, 6u};
    EXPECT_EQ(cr_axis.range(2.5, nhood11s), expected_range);
    expected_range = {0u, 9u};
    EXPECT_EQ(cr_axis.range(2.5, nhoodAlls), expected_range);
}

TEST(grid, circular_regular_axis) {

    scalar epsilon = 10 * std::numeric_limits<scalar>::epsilon();

    // Let's say 36 modules, but with 4 directly at 0, pi/2, pi, -pi2
    scalar pi{M_PI};
    scalar pi2{M_PI_2};

    scalar half_module{scalar{2} * pi2 / scalar{72}};
    scalar phi_min{-pi + half_module};
    scalar phi_max{pi - half_module};

    // Lower bin edges: min and max bin edge for the regular axis
    vecmem::vector<scalar> bin_edges = {-10, phi_min, phi_max, 56};
    // Regular axis: first entry is the offset in the bin edges vector (1), the
    // second entry is the number of bins (36)
    dindex_range edge_range = {1, 36};

    // A closed regular x-axis
    single_axis<circular<>, regular<>> cr_axis(&edge_range, &bin_edges);

    // Test axis shape
    EXPECT_EQ(cr_axis.label(), n_axis::label::e_phi);
    EXPECT_EQ(cr_axis.shape(), n_axis::shape::e_circular);
    EXPECT_EQ(cr_axis.binning(), n_axis::binning::e_regular);

    // N bins
    EXPECT_EQ(cr_axis.nbins(), 36u);
    // Axis bin access
    // overflow
    EXPECT_EQ(cr_axis.bin(phi_max + epsilon), 0u);
    // underflow
    EXPECT_EQ(cr_axis.bin(phi_min - epsilon), 35u);
    // middle of the axis
    EXPECT_EQ(cr_axis.bin(0), 19u);

    // Bin wrapping test
    typename single_axis<circular<>, regular<>>::shape_type circ_shape{};
    EXPECT_EQ(circ_shape.wrap(4, 36u), 3u);
    EXPECT_EQ(circ_shape.wrap(0, 36u), 35u);
    EXPECT_EQ(circ_shape.wrap(-1, 36u), 34u);
    EXPECT_EQ(circ_shape.wrap(36, 36u), 35u);
    EXPECT_EQ(circ_shape.wrap(40, 36u), 3u);

    // Axis range access - binned (symmetric & asymmetric)
    darray<dindex, 2> nhood00i = {0u, 0u};
    darray<dindex, 2> nhood01i = {0u, 1u};
    darray<dindex, 2> nhood11i = {1u, 1u};
    darray<dindex, 2> nhood22i = {2u, 2u};

    dindex_range expected_range = {0u, 0u};
    EXPECT_EQ(cr_axis.range(pi + epsilon, nhood00i), expected_range);
    expected_range = {0u, 1u};
    EXPECT_EQ(cr_axis.range(pi + epsilon, nhood01i), expected_range);
    expected_range = {35u, 1u};
    EXPECT_EQ(cr_axis.range(pi + epsilon, nhood11i), expected_range);
    expected_range = {34u, 2u};
    EXPECT_EQ(cr_axis.range(pi + epsilon, nhood22i), expected_range);

    // Axis range access - scalar (symmetric & asymmteric)
    darray<scalar, 2> nhood00s = {0., 0.};
    darray<scalar, 2> nhoodEpsilon = {scalar{0.5} * epsilon,
                                      scalar{0.5} * epsilon};
    scalar bin_step{cr_axis.bin_width()};
    darray<scalar, 2> nhood22s = {2 * bin_step, 2 * bin_step};

    expected_range = {0u, 0u};
    EXPECT_EQ(cr_axis.range(pi + epsilon, nhood00s), expected_range);
    EXPECT_EQ(cr_axis.range(pi + epsilon, nhoodEpsilon), expected_range);
    expected_range = {34u, 2u};
    EXPECT_EQ(cr_axis.range(pi + epsilon, nhood22s), expected_range);
}

TEST(grid, closed_irregular_axis) {

    // Lower bin edges: all lower bin edges for irregular binning, plus the
    // final upper bin edge
    vecmem::vector<scalar> bin_edges = {-100, -3, 1, 2, 4, 8, 12, 15, 18};
    // Index range for the bin edges [-3, 15]
    dindex_range edge_range = {1, 7};

    // A closed irregular z-axis
    single_axis<closed<label::e_z>, irregular<>> cir_axis(&edge_range,
                                                          &bin_edges);

    // Test axis shape
    EXPECT_EQ(cir_axis.label(), n_axis::label::e_z);
    EXPECT_EQ(cir_axis.shape(), n_axis::shape::e_closed);
    EXPECT_EQ(cir_axis.binning(), n_axis::binning::e_irregular);

    // Axis bin access
    //
    // N bins
    EXPECT_EQ(cir_axis.nbins(), 6u);
    // Bin tests
    EXPECT_EQ(cir_axis.bin(-2), 0u);
    EXPECT_EQ(cir_axis.bin(10), 4u);
    EXPECT_EQ(cir_axis.bin(5.8), 3u);
    // Underflow test
    EXPECT_EQ(cir_axis.bin(-4), 0u);
    // Overflow test
    EXPECT_EQ(cir_axis.bin(17), 5u);

    // Axis range access - binned  (symmetric & asymmetric)
    darray<dindex, 2> nhood00i = {0u, 0u};
    darray<dindex, 2> nhood01i = {0u, 1u};
    darray<dindex, 2> nhood11i = {1u, 1u};
    darray<dindex, 2> nhood22i = {2u, 2u};

    dindex_range expected_range = {2u, 2u};
    EXPECT_EQ(cir_axis.range(3., nhood00i), expected_range);
    expected_range = {1u, 3u};
    EXPECT_EQ(cir_axis.range(3., nhood11i), expected_range);
    expected_range = {2u, 3u};
    EXPECT_EQ(cir_axis.range(3., nhood01i), expected_range);
    expected_range = {0u, 1u};
    EXPECT_EQ(cir_axis.range(0., nhood11i), expected_range);
    expected_range = {2u, 5u};
    EXPECT_EQ(cir_axis.range(10., nhood22i), expected_range);

    // Axis range access - scalar
    darray<scalar, 2> nhood00s = {0., 0.};
    darray<scalar, 2> nhood10s = {1.5, 0.2};
    darray<scalar, 2> nhood11s = {4., 5.5};

    expected_range = {2u, 2u};
    EXPECT_EQ(cir_axis.range(3., nhood00s), expected_range);
    expected_range = {1u, 2u};
    EXPECT_EQ(cir_axis.range(3., nhood10s), expected_range);
    expected_range = {0u, 4u};
    EXPECT_EQ(cir_axis.range(3., nhood11s), expected_range);
}

TEST(grid, multi_axis) {

    // readable axis ownership definition
    bool constexpr is_owning = true;
    bool constexpr is_not_owning = false;

    // Lower bin edges for all axes
    vecmem::vector<scalar> bin_edges = {-10., 10., -20., 20., 0., 100.};
    // Offsets into edges container and #bins for all axes
    vecmem::vector<dindex_range> edge_ranges = {
        {0u, 20u}, {2u, 40u}, {4u, 50u}};

    // Non-owning multi axis test
    cartesian_3D<is_not_owning, host_container_types> axes(&edge_ranges,
                                                           &bin_edges);

    EXPECT_EQ(axes.Dim, 3u);

    // Get single axis objects
    auto x_axis = axes.get_axis<label::e_x>();
    EXPECT_EQ(x_axis.nbins(), 20u);
    auto y_axis = axes.get_axis<label::e_y>();
    EXPECT_EQ(y_axis.nbins(), 40u);
    auto z_axis = axes.get_axis<label::e_z>();
    EXPECT_EQ(z_axis.nbins(), 50u);

    auto ax = axes.get_axis<decltype(y_axis)>();
    EXPECT_EQ(ax.nbins(), 40u);

    // Test bin search
    point3 p3{0., 0., 0.};  // origin
    dmulti_index<dindex, 3> expected_bins{10u, 20u, 0u};
    EXPECT_EQ(axes.bins(p3), expected_bins);
    p3 = {-5.5, -3.2, 24.1};
    expected_bins = {4u, 16u, 12u};
    EXPECT_EQ(axes.bins(p3), expected_bins);
    p3 = {-1., 21., 12.};
    expected_bins = {9u, 39u, 6u};
    EXPECT_EQ(axes.bins(p3), expected_bins);

    // Test bin range search
    p3 = {1., 1., 10.};
    // Axis range access - binned  (symmetric & asymmetric)
    darray<dindex, 2> nhood00i = {0u, 0u};
    darray<dindex, 2> nhood01i = {0u, 1u};
    darray<dindex, 2> nhood22i = {2u, 2u};

    dmulti_index<dindex_range, 3> expected_ranges{};
    expected_ranges[0] = {11u, 11u};
    expected_ranges[1] = {21u, 21};
    expected_ranges[2] = {5u, 5u};
    EXPECT_EQ(axes.bin_ranges(p3, nhood00i), expected_ranges);
    expected_ranges[0] = {11u, 12u};
    expected_ranges[1] = {21u, 22u};
    expected_ranges[2] = {5u, 6u};
    EXPECT_EQ(axes.bin_ranges(p3, nhood01i), expected_ranges);
    expected_ranges[0] = {9u, 13u};
    expected_ranges[1] = {19u, 23u};
    expected_ranges[2] = {3u, 7u};
    EXPECT_EQ(axes.bin_ranges(p3, nhood22i), expected_ranges);

    // Axis range access - scalar
    darray<scalar, 2> nhood00s = {0., 0.};
    darray<scalar, 2> nhood10s = {1.5, 0.2};
    darray<scalar, 2> nhood11s = {4., 5.5};

    expected_ranges[0] = {11u, 11u};
    expected_ranges[1] = {21u, 21u};
    expected_ranges[2] = {5u, 5u};
    EXPECT_EQ(axes.bin_ranges(p3, nhood00s), expected_ranges);
    expected_ranges[0] = {9u, 11u};
    expected_ranges[1] = {19u, 21u};
    expected_ranges[2] = {4u, 5u};
    EXPECT_EQ(axes.bin_ranges(p3, nhood10s), expected_ranges);
    expected_ranges[0] = {7u, 16u};
    expected_ranges[1] = {17u, 26u};
    expected_ranges[2] = {3u, 7u};
    EXPECT_EQ(axes.bin_ranges(p3, nhood11s), expected_ranges);

    // Owning multi axis test
    vecmem::vector<scalar> bin_edges_cp(bin_edges);
    // Offsets into edges container and #bins for all axes
    vecmem::vector<dindex_range> edge_ranges_cp(edge_ranges);

    cartesian_3D<is_owning, host_container_types> axes_own(
        std::move(edge_ranges_cp), std::move(bin_edges_cp));

    EXPECT_EQ(axes_own.bins(p3), axes.bins(p3));
    EXPECT_EQ(axes_own.bin_ranges(p3, nhood00i), axes.bin_ranges(p3, nhood00i));
    EXPECT_EQ(axes_own.bin_ranges(p3, nhood01i), axes.bin_ranges(p3, nhood01i));
    EXPECT_EQ(axes_own.bin_ranges(p3, nhood22i), axes.bin_ranges(p3, nhood22i));
    EXPECT_EQ(axes_own.bin_ranges(p3, nhood00s), axes.bin_ranges(p3, nhood00s));
    EXPECT_EQ(axes_own.bin_ranges(p3, nhood10s), axes.bin_ranges(p3, nhood10s));
    EXPECT_EQ(axes_own.bin_ranges(p3, nhood11s), axes.bin_ranges(p3, nhood11s));

    // Transfer to an owning multi-axis, which uses vecmem::device_vector
    cartesian_3D<is_owning, host_container_types>::view_type axes_view =
        get_data(axes_own);
    cartesian_3D<is_owning, device_container_types> axes_device(axes_view);

    // Get single axis objects
    auto x_axis_device = axes_device.get_axis<label::e_x>();
    EXPECT_EQ(x_axis_device.nbins(), 20u);
    auto y_axis_device = axes_device.get_axis<label::e_y>();
    EXPECT_EQ(y_axis_device.nbins(), 40u);
    auto z_axis_device = axes_device.get_axis<label::e_z>();
    EXPECT_EQ(z_axis_device.nbins(), 50u);
}
