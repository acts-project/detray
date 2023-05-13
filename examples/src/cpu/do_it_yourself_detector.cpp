/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/common/detector_writer.hpp"
#include "detray/tools/surface_factory.hpp"
#include "detray/tools/volume_builder.hpp"

// Example include(s)
#include "detray/examples/detector_metadata.hpp"
#include "detray/examples/square_surface_generator.hpp"
#include "detray/examples/types.hpp"  // linear algebra types

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <memory>

/// Write a dector using the json IO
int main() {

    // The new detector type
    using detector_t = detray::detector<detray::example::my_metadata>;

    // First, create an empty detector in in host memory to be filled
    vecmem::host_memory_resource host_mr;
    detector_t det{host_mr};

    // Now fill the detector

    // Get a generic volume builder first and decorate it later
    detray::volume_builder<detector_t> vbuilder{};

    // Create a new cuboid volume in the detector
    vbuilder.init_vol(det, detray::volume_id::e_cuboid);

    // Fill some squares into the volume
    using square_factory_t =
        detray::surface_factory<detector_t, detray::square2D<>,
                                detector_t::masks::id::e_square2,
                                detray::surface_id::e_sensitive>;
    auto square_factory = std::make_shared<square_factory_t>();

    // Add a square that is 20x20mm large, likns back to its mother volume (0)
    // and is placed with a translation of (x = 1., y = 2., z = 3.)
    detray::example::vector3 translation{
        1.f * detray::unit<detray::scalar>::mm,
        2.f * detray::unit<detray::scalar>::mm,
        3.f * detray::unit<detray::scalar>::mm};
    square_factory->push_back({detray::example::transform3{translation},
                               0u,
                               {20.f * detray::unit<detray::scalar>::mm}});

    // Add some programmatically generated surfaces
    auto sq_generator =
        std::make_shared<detray::example::square_surface_generator>(
            10, 10.f * detray::unit<detray::scalar>::mm);

    // Add surfaces to volume and add the volume to the detector
    vbuilder.add_sensitives(square_factory);
    vbuilder.add_sensitives(sq_generator);

    vbuilder.build(det);

    // Write the detector to file
    auto writer_cfg = detray::io::detector_writer_config{}
                          .format(detray::io::format::json)
                          .replace_files(true);
    detray::io::write_detector(det, {{0u, "example_detector"}}, writer_cfg);
}
