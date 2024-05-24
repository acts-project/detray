/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/frontend/detector_reader_config.hpp"
#include "detray/io/frontend/detector_writer_config.hpp"
#include "detray/options/options_handling.hpp"

// Boost
#include <boost/program_options.hpp>

// System include(s)
#include <stdexcept>
#include <string>

namespace detray::options {

/// Add options for the detray detector reader
template <>
void add_options<detray::io::detector_reader_config>(
    boost::program_options::options_description &desc,
    const detray::io::detector_reader_config &) {

    desc.add_options()("geometry_file",
                       boost::program_options::value<std::string>(),
                       "Detector geometry input file")(
        "grid_file", boost::program_options::value<std::string>(),
        "Detector surface grid input file")(
        "material_file", boost::program_options::value<std::string>(),
        "Detector material input file");
}

/// Configure the detray detector reader
template <>
void configure_options<detray::io::detector_reader_config>(
    boost::program_options::variables_map &vm,
    detray::io::detector_reader_config &cfg) {

    // Input files
    if (vm.count("geometry_file")) {
        cfg.add_file(vm["geometry_file"].as<std::string>());
    } else {
        throw std::invalid_argument(
            "Please specify a geometry input file!\n\n");
    }
    if (vm.count("material_file")) {
        cfg.add_file(vm["material_file"].as<std::string>());
    }
    if (vm.count("grid_file")) {
        cfg.add_file(vm["grid_file"].as<std::string>());
    }
}

/// Add options for the detray detector writer
template <>
void add_options<detray::io::detector_writer_config>(
    boost::program_options::options_description &desc,
    const detray::io::detector_writer_config &cfg) {

    desc.add_options()(
        "outdir",
        boost::program_options::value<std::string>()->default_value(cfg.path()),
        "Output directory for detector files")("compactify_json",
                                               "not implemented")(
        "write_material", "toggle material output")("write_grids",
                                                    "toggle grid output");
}

/// Configure the detray detector writer
template <>
void configure_options<detray::io::detector_writer_config>(
    boost::program_options::variables_map &vm,
    detray::io::detector_writer_config &cfg) {

    if (!vm["outdir"].defaulted()) {
        cfg.path(vm["outdir"].as<std::string>());
    }
    if (vm.count("compactify_json")) {
        cfg.compactify_json(true);
    }
    if (vm.count("write_material")) {
        cfg.write_material(true);
    }
    if (vm.count("write_grids")) {
        cfg.write_grids(true);
    }
}

}  // namespace detray::options
