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
#include <filesystem>
#include <ios>

using namespace detray;

namespace {

/// Compare two files with names @param file_name1 and @param file_name2 for
/// equality, while skipping the first @param skip lines (header part)
bool compare_files(const std::string& file_name1, const std::string& file_name2,
                   std::size_t skip = 15u) {
    auto file1 = io::detail::file_handle(
        file_name1, std::ios_base::in | std::ios_base::binary);
    auto file2 = io::detail::file_handle(
        file_name2, std::ios_base::in | std::ios_base::binary);

    std::string line1, line2;

    // Check files line by line
    std::size_t i{1u};
    while (std::getline(*file1, line1)) {
        if (std::getline(*file2, line2)) {
            if (skip < i and line1 != line2) {
                std::cout << "In line " << i << ":" << std::endl
                          << line1 << std::endl
                          << line2 << std::endl;
                return false;
            }
        } else {
            std::cout << "Could not read next line from file 2:" << std::endl
                      << "In line " << i << ":" << std::endl
                      << line1 << std::endl;
            return false;
        }
        ++i;
    }

    // Are there more lines in file2 than file1?
    if (std::getline(*file2, line2)) {
        std::cout << "Could not read next line from file 1:" << std::endl
                  << "In line " << i << ":" << std::endl
                  << line2 << std::endl;
        return false;
    }

    // Passed
    return true;
}

}  // anonymous namespace

/// Test the reading and writing of a telescope detector
TEST(io, json_telescope_detector_reader) {

    mask<rectangle2D<>> rec2{0u, 100.f, 100.f};

    // Surface positions
    std::vector<scalar> positions = {1.f,   50.f,  100.f, 150.f, 200.f, 250.f,
                                     300.f, 350.f, 400.f, 450.f, 500.f};

    tel_det_config<rectangle2D<>> tel_cfg{rec2};
    tel_cfg.positions(positions);

    // Telescope detector
    vecmem::host_memory_resource host_mr;
    auto [telescope_det, telescope_names] =
        create_telescope_detector(host_mr, tel_cfg);

    using detector_t = decltype(telescope_det);

    auto writer_cfg = io::detector_writer_config{}
                          .format(io::format::json)
                          .write_grids(false)
                          .replace_files(true);
    io::write_detector(telescope_det, telescope_names, writer_cfg);

    // Read the detector back in
    io::detector_reader_config reader_cfg{};
    reader_cfg.do_check(true)
        .add_file("telescope_detector_geometry.json")
        .add_file("telescope_detector_homogeneous_material.json");

    const auto [det, names] =
        io::read_detector<detector_t>(host_mr, reader_cfg);

    const auto mat_store = det.material_store();
    const auto slabs = mat_store.get<detector_t::materials::id::e_slab>();

    EXPECT_EQ(det.volumes().size(), 1u);
    EXPECT_EQ(slabs.size(), positions.size() + 6u);

    // Write the result to a different set of files
    writer_cfg.replace_files(false);
    io::write_detector(det, names, writer_cfg);

    // Compare writing round-trip
    EXPECT_TRUE(compare_files("telescope_detector_geometry.json",
                              "telescope_detector_geometry_2.json"));
    EXPECT_TRUE(
        compare_files("telescope_detector_homogeneous_material.json",
                      "telescope_detector_homogeneous_material_2.json"));

    // Remove files
    std::filesystem::remove("telescope_detector_geometry_2.json");
    std::filesystem::remove("telescope_detector_homogeneous_material_2.json");
}

/// Test the reading and writing of a toy detector geometry
TEST(io, json_toy_geometry) {

    using detector_t = detector<toy_metadata>;

    // Toy detector
    vecmem::host_memory_resource host_mr;
    toy_det_config toy_cfg{};
    toy_cfg.use_material_maps(false);
    auto [toy_det, names] = create_toy_geometry(host_mr, toy_cfg);

    // Write the detector
    json_geometry_writer<detector_t> geo_writer;
    auto file_name = geo_writer.write(
        toy_det, names, std::ios::out | std::ios::binary | std::ios::trunc);

    // Empty volume name map to be filled
    typename detector_t::name_map volume_name_map = {{0u, "toy_detector"}};

    // Read the detector back in
    detector_builder<toy_metadata> toy_builder;
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
    EXPECT_EQ(comp_det.surfaces().size(), 3244u);
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

    using detector_t = detector<toy_metadata>;

    // Toy detector
    vecmem::host_memory_resource host_mr;
    toy_det_config toy_cfg{};
    toy_cfg.use_material_maps(false);
    const auto [toy_det, toy_names] = create_toy_geometry(host_mr, toy_cfg);

    auto writer_cfg = io::detector_writer_config{}
                          .format(io::format::json)
                          .replace_files(true);
    io::write_detector(toy_det, toy_names, writer_cfg);

    // Read the detector back in
    io::detector_reader_config reader_cfg{};
    reader_cfg.do_check(true)
        .add_file("toy_detector_geometry.json")
        .add_file("toy_detector_homogeneous_material.json")
        .add_file("toy_detector_surface_grids.json");

    const auto [det, names] =
        io::read_detector<detector_t, 1u>(host_mr, reader_cfg);

    EXPECT_TRUE(test_toy_detector(det, names));

    // Write the result to a different set of files
    writer_cfg.replace_files(false);
    io::write_detector(det, names, writer_cfg);

    // Compare writing round-trip
    EXPECT_TRUE(compare_files("toy_detector_geometry.json",
                              "toy_detector_geometry_2.json"));
    EXPECT_TRUE(compare_files("toy_detector_homogeneous_material.json",
                              "toy_detector_homogeneous_material_2.json"));
    EXPECT_TRUE(compare_files("toy_detector_surface_grids.json",
                              "toy_detector_surface_grids_2.json"));

    // Remove files
    std::filesystem::remove("toy_detector_geometry_2.json");
    std::filesystem::remove("toy_detector_homogeneous_material_2.json");
    std::filesystem::remove("toy_detector_surface_grids_2.json");
}

/// Test the reading and writing of a wire chamber
TEST(io, json_wire_chamber_reader) {

    using detector_t = detector<default_metadata>;

    // Wire chamber
    vecmem::host_memory_resource host_mr;
    wire_chamber_config wire_cfg{};
    wire_cfg.use_material_maps(false);
    auto [wire_det, wire_names] = create_wire_chamber(host_mr, wire_cfg);

    auto writer_cfg = io::detector_writer_config{}
                          .format(io::format::json)
                          .replace_files(true);
    io::write_detector(wire_det, wire_names, writer_cfg);

    // Read the detector back in
    io::detector_reader_config reader_cfg{};
    reader_cfg.do_check(true)
        .add_file("wire_chamber_geometry.json")
        .add_file("wire_chamber_homogeneous_material.json")
        .add_file("wire_chamber_surface_grids.json");

    const auto [det, names] =
        io::read_detector<detector_t>(host_mr, reader_cfg);

    EXPECT_EQ(det.volumes().size(), 11u);

    // Write the result to a different set of files
    writer_cfg.replace_files(false);
    io::write_detector(det, names, writer_cfg);

    // Compare writing round-trip
    EXPECT_TRUE(compare_files("wire_chamber_geometry.json",
                              "wire_chamber_geometry_2.json"));
    EXPECT_TRUE(compare_files("wire_chamber_homogeneous_material.json",
                              "wire_chamber_homogeneous_material_2.json"));
    EXPECT_TRUE(compare_files("wire_chamber_surface_grids.json",
                              "wire_chamber_surface_grids_2.json"));

    // Remove files
    std::filesystem::remove("wire_chamber_geometry_2.json");
    std::filesystem::remove("wire_chamber_homogeneous_material_2.json");
    std::filesystem::remove("wire_chamber_surface_grids_2.json");
}
