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

using namespace detray;

TEST(tools, volume_radiation_length) {

    // Memory resource
    vecmem::host_memory_resource host_mr;

    // Use rectangular surfaces
    constexpr bool rectangular = false;

    // Build from given module positions
    std::vector<scalar> positions = {0.,   50., 100., 150., 200., 250.,
                                     300., 350, 400,  450., 500.};

    // Build telescope detector with unbounded planes
    const auto det = create_telescope_detector<rectangular>(host_mr, positions);

    volume_radiation_length<decltype(det)> vrl(det);
}