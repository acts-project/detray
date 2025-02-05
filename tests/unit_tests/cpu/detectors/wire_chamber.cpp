// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

// Project include(s).
#include "detray/utils/consistency_checker.hpp"

// Detray test include(s)
#include "detray/test/utils/detectors/build_wire_chamber.hpp"
#include "detray/test/utils/types.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include.
#include <gtest/gtest.h>

using namespace detray;

GTEST_TEST(detray_detectors, wire_chamber) {

    vecmem::host_memory_resource host_mr;

    wire_chamber_config<test::scalar> cfg{};
    auto [wire_det, names] = build_wire_chamber<test::algebra>(host_mr, cfg);

    // Check general consistency of the detector
    detail::check_consistency(wire_det, true, names);
}
