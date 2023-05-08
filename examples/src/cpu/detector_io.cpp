/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/examples/types.hpp"
#include "detray/io/json/json_reader.hpp"
#include "detray/io/json/json_writer.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <ios>

/// Read and write a dector using the json IO
main() {

    // Toy detector
    using detector_t =
        detray::detector<detray::detector_registry::toy_detector>;

    vecmem::host_memory_resource host_mr;
    detector_t toy_det = detray::create_toy_geometry(host_mr);

    // Names for volumes
    typename detector_t::name_map volume_name_map = {{0u, "toy_detector"}};

    // Write the detector geometry
    detray::json_geometry_writer<detector_t> geo_writer;
    auto file_name = geo_writer.write(
        toy_det, volume_name_map, std::ios_base::out | std::ios_base::trunc);

    // Read the detector geometry back in
    detector_t det{host_mr};
    detray::json_geometry_reader<detector_t> geo_reader;
    geo_reader.read(det, volume_name_map, file_name);
}
