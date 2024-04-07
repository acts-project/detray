/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/test/detail/register_checks.hpp"
#include "detray/test/detector_consistency.hpp"
#include "detray/test/detector_helix_scan.hpp"
#include "detray/test/detector_ray_scan.hpp"
#include "detray/test/helix_navigation.hpp"
#include "detray/test/straight_line_navigation.hpp"

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
    po::options_description desc("\ndetray detector validation options");

    std::vector<dindex> window;
    desc.add_options()("help", "produce help message")(
        "write_volume_graph", "writes the volume graph to file")(
        "write_scan_data", "writes the ray/helix scan intersections to file")(
        "geometry_file", po::value<std::string>(), "geometry input file")(
        "grid_file", po::value<std::string>(), "Surface grid input file")(
        "material_file", po::value<std::string>(), "material input file")(
        "phi_steps", po::value<std::size_t>()->default_value(50u),
        "# phi steps for particle gun")(
        "theta_steps", po::value<std::size_t>()->default_value(50u),
        "# theta steps for particle gun")(
        "theta_range", po::value<std::vector<scalar_t>>()->multitoken(),
        "min, max range of theta values for particle gun")(
        "origin", po::value<std::vector<scalar_t>>()->multitoken(),
        "coordintates for particle gun origin position")(
        "p_tot", po::value<scalar_t>()->default_value(10.f),
        "total momentum of the test particle [GeV]")(
        "p_T", po::value<scalar_t>()->default_value(10.f),
        "transverse momentum of the test particle [GeV]")(
        "search_window", po::value<std::vector<dindex>>(&window)->multitoken(),
        "search window size for the grid")(
        "overstep_tol", po::value<scalar_t>()->default_value(-100.f),
        "overstepping tolerance [um] NOTE: Must be negative!");

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

    // Configs to be filled
    detray::io::detector_reader_config reader_cfg{};
    reader_cfg.do_check(false);  // < Don't run consistency check twice
    detray::test::consistency_check<detector_t>::config con_chk_cfg{};
    detray::test::ray_scan<detector_t>::config ray_scan_cfg{};
    detray::test::helix_scan<detector_t>::config hel_scan_cfg{};
    detray::test::straight_line_navigation<detector_t>::config str_nav_cfg{};
    detray::test::helix_navigation<detector_t>::config hel_nav_cfg{};

    // General options
    if (vm.count("write_volume_graph")) {
        con_chk_cfg.write_graph(true);
        throw std::invalid_argument("Writing of voume graph not implemented");
    }
    if (vm.count("write_scan_data")) {
        ray_scan_cfg.write_intersections(true);
        hel_scan_cfg.write_intersections(true);
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
    if (vm.count("grid_file")) {
        reader_cfg.add_file(vm["grid_file"].as<std::string>());
    }

    // Particle gun
    if (vm.count("phi_steps")) {
        const std::size_t phi_steps{vm["phi_steps"].as<std::size_t>()};

        ray_scan_cfg.track_generator().phi_steps(phi_steps);
        hel_scan_cfg.track_generator().phi_steps(phi_steps);
    }
    if (vm.count("theta_steps")) {
        const std::size_t theta_steps{vm["theta_steps"].as<std::size_t>()};

        ray_scan_cfg.track_generator().theta_steps(theta_steps);
        hel_scan_cfg.track_generator().theta_steps(theta_steps);
    }
    if (vm.count("theta_range")) {
        const auto theta_range = vm["theta_range"].as<std::vector<scalar_t>>();
        if (theta_range.size() == 2u) {
            ray_scan_cfg.track_generator().theta_range(theta_range[0],
                                                       theta_range[1]);
            hel_scan_cfg.track_generator().theta_range(theta_range[0],
                                                       theta_range[1]);
        } else {
            throw std::invalid_argument("Theta range needs two arguments");
        }
    }
    if (vm.count("origin")) {
        const auto origin = vm["origin"].as<std::vector<scalar_t>>();
        if (origin.size() == 3u) {
            ray_scan_cfg.track_generator().origin(
                {origin[0], origin[1], origin[2]});
            hel_scan_cfg.track_generator().origin(
                {origin[0], origin[1], origin[2]});
        } else {
            throw std::invalid_argument(
                "Particle gun origin needs three arguments");
        }
    }
    if (vm.count("p_tot")) {
        const scalar_t p_mag{vm["p_tot"].as<scalar_t>()};

        hel_scan_cfg.track_generator().p_tot(p_mag * unit<scalar_t>::GeV);
    } else if (vm.count("p_T")) {
        const scalar_t p_T{vm["p_T"].as<scalar_t>()};

        hel_scan_cfg.track_generator().p_T(p_T * unit<scalar_t>::GeV);
    }

    // Navigation
    // Grid neighborhood size
    if (vm.count("search_window")) {
        if (window.size() != 2u) {
            throw std::invalid_argument(
                "Incorrect surface grid search window. Please provide two "
                "integer distances.");
        }
        str_nav_cfg.propagation().navigation.search_window = {window[0],
                                                              window[1]};
        hel_nav_cfg.propagation().navigation.search_window = {window[0],
                                                              window[1]};
    }

    if (vm.count("overstep_tol")) {
        const scalar_t overstep_tol{vm["overstep_tol"].as<scalar_t>()};

        hel_nav_cfg.propagation().navigation.overstep_tolerance =
            overstep_tol * unit<scalar_t>::um;
    }

    vecmem::host_memory_resource host_mr;

    const auto [det, names] =
        detray::io::read_detector<detector_t>(host_mr, reader_cfg);

    // General data consistency of the detector
    detray::detail::register_checks<detray::test::consistency_check>(
        det, names, con_chk_cfg);

    // Navigation link consistency, discovered by ray intersection
    ray_scan_cfg.name("ray_scan_" + names.at(0));
    detray::detail::register_checks<detray::test::ray_scan>(det, names,
                                                            ray_scan_cfg);

    // Navigation link consistency, discovered by helix intersection
    hel_scan_cfg.name("helix_scan_" + names.at(0));
    detray::detail::register_checks<detray::test::helix_scan>(det, names,
                                                              hel_scan_cfg);

    // Comparision of straight line navigation with ray scan
    str_nav_cfg.name("straight_line_navigation_" + names.at(0));
    detray::detail::register_checks<detray::test::straight_line_navigation>(
        det, names, str_nav_cfg);

    // Comparision of navigation in a constant B-field with helix
    hel_nav_cfg.name("helix_navigation_" + names.at(0));
    detray::detail::register_checks<detray::test::helix_navigation>(
        det, names, hel_nav_cfg);

    // Run the checks
    return RUN_ALL_TESTS();
}
