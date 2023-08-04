/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/detectors/create_wire_chamber.hpp"
#include "detray/io/common/detector_reader.hpp"
#include "detray/io/common/detector_writer.hpp"
#include "detray/io/json/json_reader.hpp"
#include "detray/io/json/json_writer.hpp"
#include "detray/utils/consistency_checker.hpp"
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

    detail::check_consistency(comp_det);
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

/// Test the reading and writing of a telescope detector
TEST(io, json_telescope_detector_reader) {

    mask<rectangle2D<>> rec2{0u, 100.f, 100.f};

    // Surface positions
    std::vector<scalar> positions = {0.f,   50.f,  100.f, 150.f, 200.f, 250.f,
                                     300.f, 350.f, 400.f, 450.f, 500.f};

    tel_det_config<rectangle2D<>> tel_cfg{rec2};

    using detector_t = detector<default_metadata>;

    // Wire chamber
    vecmem::host_memory_resource host_mr;
    auto [telescope_det, telescope_names] =
        create_telescope_detector(host_mr, tel_cfg);

    auto writer_cfg = io::detector_writer_config{}
                          .format(io::format::json)
                          .replace_files(true);
    io::write_detector(telescope_det, telescope_names, writer_cfg);

    // Read the detector back in
    io::detector_reader_config reader_cfg{};
    reader_cfg.add_file("telescope_detector_geometry.json")
        .add_file("telescope_detector_homogeneous_material.json");

    const auto [det, names] =
        io::read_detector<detector_t>(host_mr, reader_cfg);

    EXPECT_EQ(det.volumes().size(), 1u);
    detray::detail::check_consistency(det);

    const auto mat_store = det.material_store();

    const auto slabs = mat_store.get<default_metadata::material_ids::e_slab>();

    EXPECT_EQ(slabs.size(), positions.size() + 6u);

    const auto surfaces = det.surface_lookup();

    for (auto& sf : surfaces) {
        std::cout << sf.material() << "  Is Sensitive: " << sf.is_sensitive()
                  << "  Thickness: " << slabs[sf.index()].thickness()
                  << std::endl;
    }
}

/// Test the reading and writing of a wire chamber
TEST(io, json_wire_chamber_reader) {

    using detector_t = detector<default_metadata>;

    // Wire chamber
    vecmem::host_memory_resource host_mr;
    auto [wire_det, wire_names] =
        create_wire_chamber(host_mr, wire_chamber_config{});

    auto writer_cfg = io::detector_writer_config{}
                          .format(io::format::json)
                          .replace_files(true);
    io::write_detector(wire_det, wire_names, writer_cfg);

    // Read the detector back in
    io::detector_reader_config reader_cfg{};
    reader_cfg.add_file("wire_chamber_geometry.json")
        .add_file("wire_chamber_homogeneous_material.json");

    const auto [det, names] =
        io::read_detector<detector_t>(host_mr, reader_cfg);

    EXPECT_EQ(det.volumes().size(), 11u);
    detray::detail::check_consistency(det);
}
