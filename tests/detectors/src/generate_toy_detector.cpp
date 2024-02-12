/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/frontend/detector_writer.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Boost
#include <boost/program_options.hpp>

int main(int argc, char **argv) {

    namespace po = boost::program_options;
    using namespace detray;

    // Options parsing
    po::options_description desc("\nToy detector generation options");

    desc.add_options()("help", "produce help message")(
        "outdir", po::value<std::string>(), "Output directory for files")(
        "write_volume_graph", "writes the volume graph to file")(
        "compactify_json", "not implemented")(
        "write_material", "toggle material output")("write_grids",
                                                    "toggle grid output")(
        "barrel_layers", po::value<unsigned int>()->default_value(4u),
        "number of barrel layers [0-4]")(
        "endcap_layers", po::value<unsigned int>()->default_value(3u),
        "number of endcap layers on either side [0-7]");

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
    detray::toy_det_config toy_cfg{};
    detray::io::detector_writer_config writer_cfg{};
    writer_cfg.format(detray::io::format::json).replace_files(false);

    // General options
    std::string outdir{vm.count("outdir") ? vm["outdir"].as<std::string>()
                                          : "./toy_detector/"};
    writer_cfg.path(std::move(outdir));
    writer_cfg.compactify_json(vm.count("compactify_json"));
    writer_cfg.write_material(vm.count("write_material"));
    writer_cfg.write_grids(vm.count("write_grids"));

    // Toy detector options
    toy_cfg.n_brl_layers(vm["barrel_layers"].as<unsigned int>());
    toy_cfg.n_edc_layers(vm["endcap_layers"].as<unsigned int>());

    // Build the geometry
    vecmem::host_memory_resource host_mr;
    auto [toy_det, toy_names] = create_toy_geometry(host_mr, toy_cfg);

    // Write to file
    detray::io::write_detector(toy_det, toy_names, writer_cfg);
}
