/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/tools/volume_radiation_length.hpp"
#include "tests/common/tools/create_telescope_detector.hpp"

// google-test include(s).
#include <gtest/gtest.h>

#include <iostream>

using namespace detray;

TEST(tools, volume_radiation_length) {

    // Memory resource
    vecmem::host_memory_resource host_mr;

    // Use rectangular surfaces
    constexpr bool rectangular = false;

    // Build from given module positions
    std::vector<scalar> positions = {0.,   50., 100., 150., 200., 250.,
                                     300., 350, 400,  450., 500.};

    // Build detector and get radiation length of detector
    const auto det_1 = create_telescope_detector<rectangular>(
        host_mr, positions, {{0, 0, 0}, 0, {0, 0, 1}, -1},
        20. * unit_constants::mm, 20. * unit_constants::mm,
        silicon_tml<scalar>(), 80 * unit_constants::um);
    volume_radiation_length<decltype(det_1)> vrl_1(det_1);

    // Same but with larger thickness
    const auto det_2 = create_telescope_detector<rectangular>(
        host_mr, positions, {{0, 0, 0}, 0, {0, 0, 1}, -1},
        20. * unit_constants::mm, 20. * unit_constants::mm,
        silicon_tml<scalar>(), 100 * unit_constants::um);
    volume_radiation_length<decltype(det_2)> vrl_2(det_2);

    EXPECT_EQ(vrl_1().size(), 1);
    EXPECT_EQ(vrl_2().size(), 1);

    // Radiation length of the first detector should be larger than the second
    // one
    EXPECT_TRUE(vrl_1[0] > vrl_2[0]);
}