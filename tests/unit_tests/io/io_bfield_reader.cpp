/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/bfield_backends.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/common/bfield_reader.hpp"
#include "detray/io/common/detector_reader.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

#include <iostream>

using namespace detray;

/// Test the reading of a magnetic field from file
TEST(io, bfield_reader) {

    // Magnetic field map using nearest neightbor interpolation
    using bfield_t = covfie::field<bfield::inhom_bknd_t>;

    // Toy detector
    using detector_t = detector<toy_metadata, bfield_t>;

    // Empty volume name map
    typename detector_t::name_map volume_name_map = {{0u, "toy_detector"}};

    // Toy detector with inhomogeneous bfield
    vecmem::host_memory_resource host_mr;
    detector_builder<toy_metadata, bfield::inhom_bknd_t> toy_builder;

    // Read the bfield from file into the detector
    covfie_reader<detector_t> bf_reader;
    bf_reader.read(toy_builder, volume_name_map,
                   !std::getenv("DETRAY_BFIELD_FILE")
                       ? ""
                       : std::getenv("DETRAY_BFIELD_FILE"));
    auto det = toy_builder.build(host_mr);

    // Test the field
    const bfield_t& bf = det.get_bfield();
    bfield_t::view_t bv(bf);

    auto field_strength = [&bv](float x, float y, float z) {
        return std::sqrt(std::pow(bv.at(x, y, z)[0], 2) +
                         std::pow(bv.at(x, y, z)[1], 2) +
                         std::pow(bv.at(x, y, z)[2], 2));
    };

    // Sample some field strengths
    constexpr scalar tol{1e-7f};
    float x = -5000.0f, y = 3200.0f, z = -7700.0f;
    ASSERT_NEAR(field_strength(x, y, z), 0.668426511f * unit<scalar>::T, tol);
    x = -2400.0f, y = 7500.0f, z = 6500.0f;
    ASSERT_NEAR(field_strength(x, y, z), 0.597449579f * unit<scalar>::T, tol);
    x = -1000.0f, y = -4300.0f, z = 2700.0f;
    ASSERT_NEAR(field_strength(x, y, z), 0.41904773f * unit<scalar>::T, tol);
    x = 6100.0f, y = 4800.0f, z = -4300.0f;
    ASSERT_NEAR(field_strength(x, y, z), 0.439638488f * unit<scalar>::T, tol);
    x = -7800.0f, y = -8400.0f, z = 9000.0f;
    ASSERT_NEAR(field_strength(x, y, z), 0.061408468f * unit<scalar>::T, tol);
    x = 6900.0f, y = -8800.0f, z = -4900.0f;
    ASSERT_NEAR(field_strength(x, y, z), 0.087212384f * unit<scalar>::T, tol);
    x = -9600.0f, y = 8600.0f, z = 3000.0f;
    ASSERT_NEAR(field_strength(x, y, z), 0.026944387f * unit<scalar>::T, tol);
    x = 6400.0f, y = -3200.0f, z = -6400.0f;
    ASSERT_NEAR(field_strength(x, y, z), 0.629775357f * unit<scalar>::T, tol);
    x = 200.0f, y = -4000.0f, z = 1700.0f;
    ASSERT_NEAR(field_strength(x, y, z), 0.004898979f * unit<scalar>::T, tol);
    x = 400.0f, y = -7100.0f, z = -11700.0f;
    ASSERT_NEAR(field_strength(x, y, z), 0.293974489f * unit<scalar>::T, tol);
}
