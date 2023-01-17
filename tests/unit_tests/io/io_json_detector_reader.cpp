/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/common/detector_reader.hpp"
#include "detray/io/common/detector_writer.hpp"
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

    // Toy detector
    vecmem::host_memory_resource host_mr;
    auto [toy_det, names] = create_toy_geometry(host_mr);

    // Write the detector
    json_geometry_writer<detector_t> geo_writer;
    auto file_name = geo_writer.write(
        toy_det, names, std::ios_base::out | std::ios_base::trunc);

    // Empty volume name map to be filled
    typename detector_t::name_map volume_name_map = {{0u, "toy_detector"}};

    // Read the detector back in
    detector_builder<toy_metadata<>, volume_builder> toy_builder;
    json_geometry_reader<detector_t> geo_reader;
    geo_reader.read(toy_builder, volume_name_map, file_name);
    auto det = toy_builder.build(host_mr);

    EXPECT_TRUE(test_toy_detector(det, volume_name_map));

    // Read the toy detector into the default detector type
    detector_builder<> comp_builder;
    json_geometry_reader<detector<>> comp_geo_reader;
    comp_geo_reader.read(comp_builder, volume_name_map, file_name);
    auto comp_det = comp_builder.build(host_mr);

    using mask_id = detector<>::masks::id;
    const auto& masks = comp_det.mask_store();

    EXPECT_EQ(comp_det.volumes().size(), 20u);
    EXPECT_EQ(comp_det.n_surfaces(), 3244u);
    EXPECT_EQ(comp_det.transform_store().size(), 3264u);
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

/// Test the reading and writing of a toy detector geometry
TEST(io, json_toy_detector_reader) {

    using detector_t = detector<toy_metadata<>>;

    // Toy detector
    vecmem::host_memory_resource host_mr;
    const auto [toy_det, toy_names] = create_toy_geometry(host_mr);

    auto writer_cfg = io::detector_writer_config{}
                          .format(io::format::json)
                          .replace_files(true);
    io::write_detector(toy_det, toy_names, writer_cfg);

    // Read the detector back in
    io::detector_reader_config reader_cfg{};
    reader_cfg.add_file("toy_detector_geometry.json")
        .add_file("toy_detector_homogeneous_material.json");

    const auto [det, names] =
        io::read_detector<detector_t>(host_mr, reader_cfg);

    EXPECT_TRUE(test_toy_detector(det, names));
}
