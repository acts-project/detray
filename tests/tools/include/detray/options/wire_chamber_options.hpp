/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray test include(s)
#include "detray/options/options_handling.hpp"
#include "detray/test/utils/detectors/create_wire_chamber.hpp"

// Boost
#include "detray/options/boost_program_options.hpp"

// System include(s)
#include <string>

namespace detray::options {

/// Add options for the detray toy detector
template <>
void add_options<wire_chamber_config>(
    boost::program_options::options_description &desc,
    const wire_chamber_config &cfg) {

    desc.add_options()(
        "layers",
        boost::program_options::value<unsigned int>()->default_value(
            cfg.n_layers()),
        "number of layers")(
        "half_z",
        boost::program_options::value<float>()->default_value(
            static_cast<float>(cfg.half_z())),
        "half length z of the chamber [mm]");
}

/// Configure the detray toy detector
template <>
void configure_options<wire_chamber_config>(
    boost::program_options::variables_map &vm, wire_chamber_config &cfg) {

    cfg.n_layers(vm["layers"].as<unsigned int>());
    cfg.half_z(vm["half_z"].as<float>());
}

}  // namespace detray::options
