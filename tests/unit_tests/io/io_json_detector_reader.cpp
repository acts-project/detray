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

    using detector_t = detector<toy_metadata<>>;

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

    // Read the toy detector into a comprehensive detector type
    detector<> comp_det{host_mr};
    json_geometry_reader<detector<>> comp_geo_reader;
    comp_geo_reader.read(comp_det, volume_name_map, file_name);

    using mask_id = detector<>::masks::id;
    const auto& masks = comp_det.mask_store();

    EXPECT_EQ(comp_det.volumes().size(), 20u);
    EXPECT_EQ(comp_det.n_surfaces(), 3244u);
    EXPECT_EQ(comp_det.transform_store().size(), 3244u);
    EXPECT_EQ(masks.template size<mask_id::e_rectangle2>(), 2492u);
    EXPECT_EQ(masks.template size<mask_id::e_portal_rectangle2>(), 2492u);
    EXPECT_EQ(masks.template size<mask_id::e_trapezoid2>(), 648u);
    EXPECT_EQ(masks.template size<mask_id::e_annulus2>(), 0u);
    EXPECT_EQ(masks.template size<mask_id::e_cylinder2>(), 0u);
    EXPECT_EQ(masks.template size<mask_id::e_portal_cylinder2>(), 52u);
    EXPECT_EQ(masks.template size<mask_id::e_ring2>(), 52u);
    EXPECT_EQ(masks.template size<mask_id::e_portal_ring2>(), 52u);
    EXPECT_EQ(masks.template size<mask_id::e_straw_wire>(), 0u);
    EXPECT_EQ(masks.template size<mask_id::e_cell_wire>(), 0u);
}
