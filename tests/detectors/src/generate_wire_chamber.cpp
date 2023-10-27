/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_wire_chamber.hpp"
#include "detray/io/common/detector_writer.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Boost
#include <boost/program_options.hpp>

int main(int argc, char **argv) {

    namespace po = boost::program_options;
    using namespace detray;

    using detector_t = detector<>;
    using scalar_t = typename detector_t::scalar_type;

    // Options parsing
    po::options_description desc("\nWire chamber generation options");

    desc.add_options()("help", "produce help message")(
        "outdir", po::value<std::string>(), "Output directory for files")(
        "write_volume_graph", "writes the volume graph to file")(
        "compactify_json", "not implemented")(
        "write_material", "toggle material output")("write_grids",
                                                    "toggle grid output")(
        "layers", po::value<unsigned int>()->default_value(10u),
        "number of layers")(
        "half_z",
        po::value<scalar_t>()->default_value(1000.f * unit<scalar_t>::mm),
        "half length z of the chamber [mm]");

    po::variables_map vm;
    po::store(parse_command_line(argc, argv, desc,
                                 po::command_line_style::unix_style ^
                                     po::command_line_style::allow_short),
              vm);
    po::notify(vm);

    // Help message
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return EXIT_FAILURE;
    }

    // Configuration
    detray::wire_chamber_config wire_cfg{};
    detray::io::detector_writer_config writer_cfg{};
    writer_cfg.format(detray::io::format::json).replace_files(false);

    // General options
    std::string outdir{vm.count("outdir") ? vm["outdir"].as<std::string>()
                                          : "./wire_chamber/"};
    writer_cfg.path(std::move(outdir));
    writer_cfg.compactify_json(vm.count("compactify_json"));
    writer_cfg.write_material(vm.count("write_material"));
    writer_cfg.write_grids(vm.count("write_grids"));

    // Wire chamber options
    wire_cfg.n_layers(vm["layers"].as<unsigned int>());
    wire_cfg.half_z(vm["half_z"].as<scalar_t>());

    // Build the geometry
    vecmem::host_memory_resource host_mr;
    auto [wire_chamber, names] = create_wire_chamber(host_mr, wire_cfg);

    // Write to file
    detray::io::write_detector(wire_chamber, names, writer_cfg);
}
