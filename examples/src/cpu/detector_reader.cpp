/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/examples/types.hpp"
#include "detray/geometry/volume_graph.hpp"
#include "detray/io/json/json_reader.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>
#include <stdexcept>
#include <string>

/// Read a detector from file. For now: Read in the geometry by calling the
/// json geometry reader directly.
int main(int argc, char** argv) {

    // Input data file
    std::string file_name;
    if (argc == 2) {
        file_name = argv[1];
    } else {
        throw std::runtime_error("Please specify an input file name!");
    }

    // Read a toy detector
    using detector_t =
        detray::detector<detray::detector_registry::toy_detector>;

    // Empty volume name map (won't be filled by the reader, yet)
    typename detector_t::name_map volume_name_map = {{0u, "toy_detector"}};

    // Create an empty detector to be filled
    vecmem::host_memory_resource host_mr;
    detector_t det{host_mr};
    // Read the json geometry file
    detray::json_geometry_reader<detector_t> geo_reader;
    geo_reader.read(det, volume_name_map, file_name);

    // The detector now contains a geometry!
    detray::volume_graph graph(det);
    std::cout << "\nRead " << det.volumes().size() << " volumes from file "
              << file_name << ":\n\n"
              << graph.to_string() << std::endl;
}
