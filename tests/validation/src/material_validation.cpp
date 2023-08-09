/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/io/common/detector_reader.hpp"
#include "detray/validation/detail/register_checks.hpp"
#include "detray/validation/detector_material_scan.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// Boost
#include <boost/program_options.hpp>

// System include(s)
#include <sstream>
#include <stdexcept>
#include <string>

namespace po = boost::program_options;
using namespace detray;

int main(int argc, char **argv) {

    // Use the most general type to be able to read in all detector files
    using detector_t = detray::detector<>;
    using scalar_t = detector_t::scalar_type;

    // Filter out the google test flags
    ::testing::InitGoogleTest(&argc, argv);

    // Options parsing
    po::options_description desc("\ndetray validation options");

    desc.add_options()("help", "produce help message")(
        "geometry_file", po::value<std::string>(), "geometry input file")(
        "material_file", po::value<std::string>(), "material input file")(
        "bfield_file", po::value<std::string>(),
        "magnetic field map input file")(
        "phi_steps", po::value<std::size_t>()->default_value(50u),
        "# phi steps for particle gun")(
        "eta_steps", po::value<std::size_t>()->default_value(50u),
        "# eta steps for particle gun");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // Help message
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return EXIT_FAILURE;
    }

    // Configs to be filled
    detray::io::detector_reader_config reader_cfg{};
    detray::material_scan<detector_t>::config mat_scan_cfg{};
    mat_scan_cfg.track_generator().uniform_eta(true);

    // Input files
    if (vm.count("geometry_file")) {
        reader_cfg.add_file(vm["geometry_file"].as<std::string>());
    } else {
        std::stringstream err_stream{};
        err_stream << "Please specify a geometry input file!\n\n" << desc;

        throw std::invalid_argument(err_stream.str());
    }
    if (vm.count("material_file")) {
        reader_cfg.add_file(vm["material_file"].as<std::string>());
    } else {
        std::stringstream err_stream{};
        err_stream << "Please specify a material input file!\n\n" << desc;

        throw std::invalid_argument(err_stream.str());
    }
    if (vm.count("bfield_file")) {
        // reader_cfg.add_file(vm["bfield_file"].as<std::string>());
        throw std::invalid_argument("Reading of bfield file not implemented");
    } else {
        reader_cfg.bfield_vec(0.f, 0.f, 2.f * unit<scalar_t>::T);
    }

    // Particle gun
    if (vm.count("phi_steps")) {
        const std::size_t phi_steps{vm["phi_steps"].as<std::size_t>()};

        mat_scan_cfg.track_generator().phi_steps(phi_steps);
    }
    if (vm.count("eta_steps")) {
        const std::size_t eta_steps{vm["eta_steps"].as<std::size_t>()};

        mat_scan_cfg.track_generator().eta_steps(eta_steps);
        mat_scan_cfg.track_generator().eta_range(-4.f, 4.f);
    }

    vecmem::host_memory_resource host_mr;

    const auto [det, names] =
        detray::io::read_detector<detector_t>(host_mr, reader_cfg);

    // Print the detector's material as recorded by a ray scan
    detray::detail::register_checks<detray::material_scan>(det, names,
                                                           mat_scan_cfg);

    // Run the checks
    return RUN_ALL_TESTS();
}
