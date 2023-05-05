/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <map>
#include <string>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/detector.hpp"
#include "detray/detectors/detector_metadata.hpp"
#include "detray/test/types.hpp"
#include "tests/common/tools/read_geometry.hpp"

/// @note test has to be defined with a preprocessor command

// This tests the construction of a detector class
GTEST_TEST(detray_core, read_detector) {
    vecmem::host_memory_resource host_mr;
    using namespace detray;

    /*auto [d, name_map] =
        read_from_csv<detector_registry::tml_detector>(tml_files, host_mr);

    std::cout << d.to_string(name_map) << std::endl;*/
}
