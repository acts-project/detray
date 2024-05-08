/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/detectors/create_wire_chamber.hpp"
#include "detray/utils/consistency_checker.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include.
#include <gtest/gtest.h>

using namespace detray;

GTEST_TEST(detray_detectors, wire_chamber) {

    vecmem::host_memory_resource host_mr;

    auto [wire_det, names] =
        create_wire_chamber(host_mr, wire_chamber_config{});

    // Check general consistency of the detector
    detail::check_consistency(wire_det, true, names);
}
