/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/test/types.hpp"
#include "tests/common/test_toy_detector.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

// This test check the building of the tml based toy geometry
GTEST_TEST(detray_detectors, toy_geometry) {

    vecmem::host_memory_resource host_mr;
    constexpr std::size_t n_brl_layers{4u};
    constexpr std::size_t n_edc_layers{3u};

    const auto [toy_det, names] =
        create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    EXPECT_TRUE(test_toy_detector(toy_det));
}
