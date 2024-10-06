/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"

// Detray IO include(s)
#include "detray/io/frontend/detector_writer.hpp"

// Detray test include(s)
#include "detray/options/detector_io_options.hpp"
#include "detray/options/parse_options.hpp"
#include "detray/options/wire_chamber_options.hpp"
#include "detray/test/utils/detectors/build_wire_chamber.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Boost
#include "detray/options/boost_program_options.hpp"

namespace po = boost::program_options;

using namespace detray;

int main(int argc, char **argv) {

    // Options parsing
    po::options_description desc("\nWire chamber generation options");

    desc.add_options()("write_volume_graph", "writes the volume graph to file");

    // Configuration
    detray::wire_chamber_config wire_cfg{};
    detray::io::detector_writer_config writer_cfg{};
    writer_cfg.format(detray::io::format::json).replace_files(false);
    // Default output path
    writer_cfg.path("./wire_chamber/");

    po::variables_map vm =
        detray::options::parse_options(desc, argc, argv, wire_cfg, writer_cfg);

    // Build the geometry
    vecmem::host_memory_resource host_mr;
    auto [wire_chamber, names] = build_wire_chamber(host_mr, wire_cfg);

    // Write to file
    detray::io::write_detector(wire_chamber, names, writer_cfg);

    // General options
    if (vm.count("write_volume_graph")) {
        throw std::invalid_argument("Writing of volume graph not implemented");
    }
}
