/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "grids/axis.hpp"
#include "tests/common/test_defs.hpp"
#include "utils/indexing.hpp"

#include <gtest/gtest.h>

#include <climits>

using namespace detray;

TEST(grids, regular_closed_axis) {

  axis::regular<> ten_bins{10, -3., 7.};
  // N bins
  EXPECT_EQ(ten_bins.bins(), 10u);
  // Axis bin access
  EXPECT_EQ(ten_bins.bin(-4.), 0u);
  EXPECT_EQ(ten_bins.bin(2.5), 5u);
  EXPECT_EQ(ten_bins.bin(8.), 9u);

  // Axis range access - binned (symmetric & asymmetric)
  darray<dindex, 2> zone00 = {0u, 0u};
  darray<dindex, 2> zone01 = {0u, 1u};
  darray<dindex, 2> zone11 = {1u, 1u};
  darray<dindex, 2> zone44 = {4u, 4u};
  darray<dindex, 2> zone55 = {5u, 5u};

  dindex_range expected_range = {5u, 5u};
  EXPECT_EQ(ten_bins.range(2.5, zone00), expected_range);
  expected_range = {4u, 6u};
  EXPECT_EQ(ten_bins.range(2.5, zone11), expected_range);
  expected_range = {5u, 6u};
  EXPECT_EQ(ten_bins.range(2.5, zone01), expected_range);
  expected_range = {0u, 8u};
  EXPECT_EQ(ten_bins.range(1.5, zone44), expected_range);
  expected_range = {3u, 9u};
  EXPECT_EQ(ten_bins.range(5.5, zone55), expected_range);

  // Axis sequence access - binned (symmetric & asymmetric)
  dindex_sequence expected_zone = {5u};
  EXPECT_EQ(ten_bins.zone(2.5, zone00), expected_zone);
  expected_zone = {5u, 6u};
  EXPECT_EQ(ten_bins.zone(2.5, zone01), expected_zone);
  expected_zone = {4u, 5u, 6u};
  EXPECT_EQ(ten_bins.zone(2.5, zone11), expected_zone);
  expected_zone = {0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u};
  EXPECT_EQ(ten_bins.zone(1.5, zone44), expected_zone);

  // Axis range access - scalar (symmteric & asymmetric)
  darray<scalar, 2> szone00 = {0., 0.};
  darray<scalar, 2> sepsilon = {0.01, 0.01};
  darray<scalar, 2> szone11 = {1., 1.};
  darray<scalar, 2> szoneAll = {10., 10.};

  expected_range = {5u, 5u};
  EXPECT_EQ(ten_bins.range(2.5, szone00), expected_range);
  EXPECT_EQ(ten_bins.range(2.5, sepsilon), expected_range);
  expected_range = {4u, 6u};
  EXPECT_EQ(ten_bins.range(2.5, szone11), expected_range);
  expected_range = {0u, 9u};
  EXPECT_EQ(ten_bins.range(2.5, szoneAll), expected_range);

  // Axis sequence acces - scalar (symmteric & asymmetric)
  expected_zone = {5u};
  EXPECT_EQ(ten_bins.zone(2.5, szone00), expected_zone);
  EXPECT_EQ(ten_bins.zone(2.5, sepsilon), expected_zone);
  expected_zone = {4u, 5u, 6u};
  EXPECT_EQ(ten_bins.zone(2.5, szone11), expected_zone);
  expected_zone = {0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u};
  EXPECT_EQ(ten_bins.zone(2.5, szoneAll), expected_zone);
}

TEST(grids, regular_circular_axis) {
  scalar epsilon = 10 * std::numeric_limits<scalar>::epsilon();

  // Let's say 36 modules, but with 4 directly at 0, pi/2, pi, -pi2
  scalar half_module = 2 * M_PI_2 / 72;
  scalar phi_min = -M_PI + half_module;
  scalar phi_max = M_PI - half_module;
  axis::circular<> full_pi = {36, phi_min, phi_max};
  // N bins
  EXPECT_EQ(full_pi.bins(), 36u);
  // Axis bin access
  EXPECT_EQ(full_pi.bin(M_PI - epsilon), 0u);
  EXPECT_EQ(full_pi.bin(M_PI + epsilon), 0u);
  EXPECT_EQ(full_pi.bin(0), 18u);
  // Remap test
  EXPECT_EQ(full_pi.remap(4, -1), 3u);
  EXPECT_EQ(full_pi.remap(4, 1), 5u);
  EXPECT_EQ(full_pi.remap(0, -1), 35);
  EXPECT_EQ(full_pi.remap(0, -2), 34);
  EXPECT_EQ(full_pi.remap(-1, -1), 34);
  EXPECT_EQ(full_pi.remap(35, 1), 0);

  // Axis range access - binned (symmetric & asymmetric)
  darray<dindex, 2> zone00 = {0u, 0u};
  darray<dindex, 2> zone01 = {0u, 1u};
  darray<dindex, 2> zone11 = {1u, 1u};
  darray<dindex, 2> zone22 = {2u, 2u};

  dindex_range expected_range = {0u, 0u};
  EXPECT_EQ(full_pi.range(M_PI + epsilon, zone00), expected_range);
  expected_range = {0u, 1u};
  EXPECT_EQ(full_pi.range(M_PI + epsilon, zone01), expected_range);
  expected_range = {35u, 1u};
  EXPECT_EQ(full_pi.range(M_PI + epsilon, zone11), expected_range);
  expected_range = {34u, 2u};
  EXPECT_EQ(full_pi.range(M_PI + epsilon, zone22), expected_range);

  // Zone test - binned
  dindex_sequence expected_zone = {34u, 35u, 0u, 1u, 2u};
  EXPECT_EQ(full_pi.zone(M_PI + epsilon, zone22), expected_zone);

  // Axis range access - scalar (symmetric & asymmteric)
  darray<scalar, 2> szone00 = {0., 0.};
  darray<scalar, 2> szoneEpsilon = {0.5 * epsilon, 0.5 * epsilon};
  scalar bin_step = (full_pi.max - full_pi.min) / full_pi.bins();
  darray<scalar, 2> szone22 = {2 * bin_step, 2 * bin_step};

  expected_range = {0u, 0u};
  EXPECT_EQ(full_pi.range(M_PI + epsilon, szone00), expected_range);
  EXPECT_EQ(full_pi.range(M_PI + epsilon, szoneEpsilon), expected_range);

  expected_range = {34u, 2u};
  EXPECT_EQ(full_pi.range(M_PI + epsilon, szone22), expected_range);

  expected_zone = {34u, 35u, 0u, 1u, 2u};
  EXPECT_EQ(full_pi.zone(M_PI + epsilon, szone22), expected_zone);
}

TEST(grids, irregular_closed_axis) {
  axis::irregular<> nonreg{{-3., 1., 2, 4., 8., 12.}};

  // Axis bin access
  //
  // N bins
  EXPECT_EQ(nonreg.bins(), 5u);
  // Bin tests
  EXPECT_EQ(nonreg.bin(-2), 0u);
  EXPECT_EQ(nonreg.bin(10), 4u);
  // Underflow test
  EXPECT_EQ(nonreg.bin(-4), 0u);
  // Overflow test
  EXPECT_EQ(nonreg.bin(14), 4u);

  // Axis range access - binned  (symmetric & asymmetric)
  darray<dindex, 2> zone00 = {0u, 0u};
  darray<dindex, 2> zone01 = {0u, 1u};
  darray<dindex, 2> zone11 = {1u, 1u};
  darray<dindex, 2> zone22 = {2u, 2u};

  dindex_range expected_range = {1u, 3u};
  EXPECT_EQ(nonreg.range(3., zone11), expected_range);

  expected_range = {2u, 3u};
  EXPECT_EQ(nonreg.range(3., zone01), expected_range);

  dindex_range expected_range_truncated_low = {0u, 1u};
  EXPECT_EQ(nonreg.range(0., zone11), expected_range_truncated_low);

  dindex_range expected_range_truncated_high = {2u, 4u};
  EXPECT_EQ(nonreg.range(10., zone22), expected_range_truncated_high);

  // Axis sequence access - binned
  dindex_sequence expected_zone = {1u, 2u, 3u};
  EXPECT_EQ(nonreg.zone(3., zone11), expected_zone);

  dindex_sequence expected_zone_truncated_low = {0u, 1u};
  EXPECT_EQ(nonreg.zone(0., zone11), expected_zone_truncated_low);

  dindex_sequence expected_zone_truncated_high = {2u, 3u, 4u};
  EXPECT_EQ(nonreg.zone(10., zone22), expected_zone_truncated_high);

  // Axis range access - scalar
  darray<scalar, 2> szone00 = {0., 0.};
  darray<scalar, 2> szone10 = {1.5, 0.2};
  expected_range = {2u, 2u};
  EXPECT_EQ(nonreg.range(3., szone00), expected_range);

  expected_range = {1u, 2u};
  EXPECT_EQ(nonreg.range(3., szone10), expected_range);

  // Axis sequence access - scalar
  expected_zone = {2u};
  EXPECT_EQ(nonreg.zone(3., szone00), expected_zone);

  expected_zone = {1u, 2u};
  EXPECT_EQ(nonreg.zone(3., szone10), expected_zone);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
