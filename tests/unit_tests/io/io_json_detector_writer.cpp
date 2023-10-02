/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/detectors/create_wire_chamber.hpp"
#include "detray/io/common/detector_writer.hpp"
#include "detray/io/json/json_writer.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <ios>

using namespace detray;

namespace {

// Use stereo annulus surfaces
constexpr scalar minR{7.2f * unit<scalar>::mm};
constexpr scalar maxR{12.0f * unit<scalar>::mm};
constexpr scalar minPhi{0.74195f};
constexpr scalar maxPhi{1.33970f};
constexpr scalar offset_x{-2.f * unit<scalar>::mm};
constexpr scalar offset_y{-2.f * unit<scalar>::mm};
mask<annulus2D<>> ann2{0u, minR, maxR, minPhi, maxPhi, offset_x, offset_y, 0.f};

// Surface positions
std::vector<scalar> positions = {0.f,   50.f,  100.f, 150.f, 200.f, 250.f,
                                 300.f, 350.f, 400.f, 450.f, 500.f};

tel_det_config<annulus2D<>> tel_cfg{ann2};

}  // anonymous namespace

/// Test the writing of a telescope detector geometry to json
TEST(io, json_telescope_geometry_writer) {

    using detector_t = detector<telescope_metadata<annulus2D<>>>;

    // Telescope detector
    vecmem::host_memory_resource host_mr;
    auto [det, names] =
        create_telescope_detector(host_mr, tel_cfg.positions(positions));

    json_geometry_writer<detector_t> geo_writer;
    geo_writer.write(det, names, std::ios::out | std::ios::binary);
}

/// Test the writing of the toy detector material to json
TEST(io, json_telescope_material_writer) {

    using detector_t = detector<telescope_metadata<annulus2D<>>>;

    // Telescope detector
    vecmem::host_memory_resource host_mr;

    auto [det, names] =
        create_telescope_detector(host_mr, tel_cfg.positions(positions));

    json_homogeneous_material_writer<detector_t> mat_writer;
    mat_writer.write(det, names, std::ios::out | std::ios::binary);
}

/// Test the writing of the toy detector grids to json
TEST(io, json_toy_grid_writer) {

    using detector_t = detector<toy_metadata>;

    // Toy detector
    vecmem::host_memory_resource host_mr;
    auto [det, names] = create_toy_geometry(host_mr);

    json_grid_writer<detector_t> grid_writer;
    grid_writer.write(det, names,
                      std::ios::out | std::ios::binary | std::ios::trunc);
}

/// Test the writing of the entire toy detector to json
TEST(io, json_toy_detector_writer) {

    // Toy detector
    vecmem::host_memory_resource host_mr;
    auto [det, names] = create_toy_geometry(host_mr);

    auto writer_cfg = io::detector_writer_config{}
                          .format(io::format::json)
                          .replace_files(true);
    io::write_detector(det, names, writer_cfg);
}

/// Test the writing of the entire wire chamber to json
TEST(io, json_wire_chamber_writer) {

    // Wire chamber
    vecmem::host_memory_resource host_mr;
    auto [det, names] = create_wire_chamber(host_mr, wire_chamber_config{});

    auto writer_cfg = io::detector_writer_config{}
                          .format(io::format::json)
                          .replace_files(true);
    io::write_detector(det, names, writer_cfg);
}
