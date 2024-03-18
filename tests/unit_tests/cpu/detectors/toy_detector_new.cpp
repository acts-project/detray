/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/detectors/build_toy_detector.hpp"
#include "detray/test/toy_detector_test_new.hpp"
#include "detray/test/types.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

// This test check the building of the tml based toy geometry
GTEST_TEST(detray_detectors, toy_detector_new) {

    vecmem::host_memory_resource host_mr;

    toy_config<scalar> toy_cfg{};
    toy_cfg.use_material_maps(false).do_check(true);
    const auto [toy_det, names] = build_toy_detector(host_mr, toy_cfg);

    EXPECT_TRUE(toy_detector_test_new(toy_det, names));

    /*toy_cfg.use_material_maps(true);
    const auto [toy_det2, names2] = build_toy_detector(host_mr, toy_cfg);

    EXPECT_TRUE(toy_detector_test_new(toy_det2, names2));*/
}
