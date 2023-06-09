/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/common/bfield_reader.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/nearest_neighbour.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>
#include <covfie/core/vector.hpp>

// GTest include(s)
#include <gtest/gtest.h>

#include <iostream>

using namespace detray;

/// Test the reading of a magnetic field from file
TEST(io, bfield_reader) {

    // Magnetic field map using nearest neightbor interpolation
    using bfield_bknd_t = covfie::backend::affine<
        covfie::backend::nearest_neighbour<covfie::backend::strided<
            covfie::vector::ulong3,
            covfie::backend::array<covfie::vector::vector_d<scalar, 3>>>>>;
    //using bfield_t = covfie::field<bfield_bknd_t>;

    // Toy detector
    using detector_t =
        detector<detector_registry::template toy_detector<bfield_bknd_t>,
                 covfie::field>;

    // @todo : Create volume name map in 'create_toy_detector'
    typename detector_t::name_map volume_name_map = {{0u, "toy_detector"}};

    // Toy detector
    vecmem::host_memory_resource host_mr;
    // Build detector and magnetic field
    detector_t toy_det = create_toy_geometry<bfield_bknd_t>(host_mr);

    // Read the bfield from file
    covfie_reader<detector_t> bf_reader;
    bf_reader.read(toy_det, volume_name_map, std::getenv("DETRAY_BFIELD_FILE"));

    // Test the field
    //const bfield_t& bf = toy_det.get_bfield();
    //bfield_t::view_t bv(bf);

    /*for (float x = -10.f; x <= 10.f; x += 1.f) {
        for (float y = -10.f; y <= 10.f; y += 1.f) {
            for (float z = -10.f; z <= 10.f; z += 1.f) {
                std::cout << bv.at(x, y, z)[0] << ", ";
                std::cout << bv.at(x, y, z)[1] << ", ";
                std::cout << bv.at(x, y, z)[2] << std::endl;
            }
        }
    }*/
}
