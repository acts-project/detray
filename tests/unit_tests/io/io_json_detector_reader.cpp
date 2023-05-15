/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/json/json_reader.hpp"
#include "detray/io/json/json_writer.hpp"
#include "tests/common/test_toy_detector.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <ios>

using namespace detray;

/// Test the reading and writing of a toy detector geometry
TEST(io, json_toy_geometry) {

    using detector_t = detector<detector_registry::toy_detector>;

    // @todo : Create volume name map in 'create_toy_detector'
    typename detector_t::name_map volume_name_map = {{0u, "toy_detector"}};

    // Toy detector
    vecmem::host_memory_resource host_mr;
    detector_t toy_det = create_toy_geometry(host_mr);

    // Write the detector
    json_geometry_writer<detector_t> geo_writer;
    auto file_name = geo_writer.write(
        toy_det, volume_name_map, std::ios_base::out | std::ios_base::trunc);

    // Read the detector back in
    detector_t det{host_mr};
    json_geometry_reader<detector_t> geo_reader;
    geo_reader.read(det, volume_name_map, file_name);

    EXPECT_TRUE(test_toy_detector(det));
}
