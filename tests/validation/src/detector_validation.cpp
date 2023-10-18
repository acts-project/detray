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
#include "detray/validation/detector_consistency.hpp"
#include "detray/validation/detector_helix_scan.hpp"
#include "detray/validation/detector_ray_scan.hpp"
#include "detray/validation/helix_navigation.hpp"
#include "detray/validation/straight_line_navigation.hpp"

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
        "write_volume_graph", "writes the volume graph to file")(
        "write_ray_scan", "writes the ray scan intersections to file")(
        "geometry_file", po::value<std::string>(), "geometry input file")(
        "material_file", po::value<std::string>(), "material input file")(
        "phi_steps", po::value<std::size_t>()->default_value(50u),
        "# phi steps for particle gun")(
        "theta_steps", po::value<std::size_t>()->default_value(50u),
        "# theta steps for particle gun")(
        "p_mag", po::value<scalar_t>()->default_value(10.f),
        "absolute momentum of the test particle [GeV]")(
        "overstep_tol", po::value<scalar_t>()->default_value(-100.f),
        "overstepping tolerance [um] NOTE: Must be negative!");

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
    detray::consistency_check<detector_t>::config con_chk_cfg{};
    detray::ray_scan<detector_t>::config ray_scan_cfg{};
    detray::helix_scan<detector_t>::config hel_scan_cfg{};
    detray::straight_line_navigation<detector_t>::config str_nav_cfg{};
    detray::helix_navigation<detector_t>::config hel_nav_cfg{};

    // General options
    if (vm.count("write_volume_graph")) {
        con_chk_cfg.write_graph(true);
        throw std::invalid_argument("Writing of voume graph not implemented");
    }
    if (vm.count("write_ray_scan")) {
        ray_scan_cfg.write_intersections(true);
        hel_scan_cfg.write_intersections(true);
        throw std::invalid_argument("Writing of ray scan not implemented");
    }

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
    }

    // Particle gun
    if (vm.count("phi_steps")) {
        const std::size_t phi_steps{vm["phi_steps"].as<std::size_t>()};

        ray_scan_cfg.track_generator().phi_steps(phi_steps);
    }
    if (vm.count("theta_steps")) {
        const std::size_t theta_steps{vm["theta_steps"].as<std::size_t>()};

        ray_scan_cfg.track_generator().theta_steps(theta_steps);
    }
    if (vm.count("p_mag")) {
        const scalar_t p_mag{vm["p_mag"].as<scalar_t>()};

        hel_nav_cfg.track_generator().p_mag(p_mag * unit<scalar_t>::GeV);
    }

    // Navigation
    if (vm.count("overstep_tol")) {
        const scalar_t overstep_tol{vm["overstep_tol"].as<scalar_t>()};

        hel_nav_cfg.overstepping_tolerance(overstep_tol * unit<scalar_t>::um);
    }

    // Ray scan and straight line navigation check use the same generator type
    str_nav_cfg.track_generator() = ray_scan_cfg.track_generator();
    hel_scan_cfg.track_generator() = hel_nav_cfg.track_generator();
    hel_scan_cfg.overstepping_tolerance(hel_nav_cfg.overstepping_tolerance());

    vecmem::host_memory_resource host_mr;

    const auto [det, names] =
        detray::io::read_detector<detector_t>(host_mr, reader_cfg);

    // General data consistency of the detector
    detray::detail::register_checks<detray::consistency_check>(det, names,
                                                               con_chk_cfg);

    // Navigation link consistency, discovered by ray intersection
    detray::detail::register_checks<detray::ray_scan>(det, names, ray_scan_cfg);

    // Navigation link consistency, discovered by helix intersection
    // Full number of helices only works in double precision
    // hel_nav_cfg.track_generator().phi_steps(100u).theta_steps(100u);
    hel_scan_cfg.track_generator().phi_steps(30u).theta_steps(30u);
    detray::detail::register_checks<detray::helix_scan>(det, names,
                                                        hel_scan_cfg);

    // Comparision of straight line navigation with ray scan
    detray::detail::register_checks<detray::straight_line_navigation>(
        det, names, str_nav_cfg);

    // Comparision of navigation in a constant B-field with helix
    // For now, run with reduced number of helices, until helix test is fixed
    hel_nav_cfg.track_generator().phi_steps(30u).theta_steps(30u);
    detray::detail::register_checks<detray::helix_navigation>(det, names,
                                                              hel_nav_cfg);

    // Run the checks
    return RUN_ALL_TESTS();
}
