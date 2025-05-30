/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"

// Detray IO include(s)
#include "detray/io/frontend/detector_reader.hpp"

// Detray test include(s)
#include "detray/options/detector_io_options.hpp"
#include "detray/options/parse_options.hpp"
#include "detray/options/propagation_options.hpp"
#include "detray/options/track_generator_options.hpp"
#include "detray/test/cpu/detector_consistency.hpp"
#include "detray/test/cpu/detector_scan.hpp"
#include "detray/test/cpu/navigation_validation.hpp"
#include "detray/test/framework/register_checks.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/framework/whiteboard.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// Boost
#include "detray/options/boost_program_options.hpp"

// System include(s)
#include <sstream>
#include <stdexcept>
#include <string>

namespace po = boost::program_options;

using namespace detray;

int main(int argc, char** argv) {

    // Use the most general type to be able to read in all detector files
    using metadata_t = test::default_metadata;
    using detector_t = detector<metadata_t>;
    using scalar = dscalar<typename detector_t::algebra_type>;

    // Filter out the google test flags
    ::testing::InitGoogleTest(&argc, argv);

    // Specific options for this test
    po::options_description desc("\ndetray detector validation options");

    desc.add_options()("write_volume_graph", "Write the volume graph to file")(
        "write_scan_data", "Write the ray/helix scan data to file")(
        "data_dir",
        po::value<std::string>()->default_value("./validation_data"),
        "Directory that contains the data files")(
        "overlaps_tol",
        po::value<float>()->default_value(stepping::config{}.min_stepsize),
        "Tolerance for considering surfaces to be overlapping [mm]");

    // Configs to be filled
    detray::io::detector_reader_config reader_cfg{};
    reader_cfg.do_check(false);  // < Don't run consistency check twice
    detray::test::consistency_check<detector_t>::config con_chk_cfg{};
    detray::test::ray_scan<detector_t>::config ray_scan_cfg{};
    detray::test::helix_scan<detector_t>::config hel_scan_cfg{};
    detray::test::straight_line_navigation<detector_t>::config str_nav_cfg{};
    detray::test::helix_navigation<detector_t>::config hel_nav_cfg{};

    po::variables_map vm = detray::options::parse_options(
        desc, argc, argv, reader_cfg, hel_scan_cfg.track_generator(),
        hel_nav_cfg.propagation());

    // General options
    if (vm.count("write_volume_graph")) {
        con_chk_cfg.write_graph(true);
        throw std::invalid_argument("Writing of volume graph not implemented");
    }
    if (vm.count("write_scan_data")) {
        ray_scan_cfg.write_intersections(true);
        hel_scan_cfg.write_intersections(true);
    }
    if (vm.count("overlaps_tol")) {
        ray_scan_cfg.overlaps_tol(vm["overlaps_tol"].as<float>());
        hel_scan_cfg.overlaps_tol(vm["overlaps_tol"].as<float>());
    }
    const auto data_dir{vm["data_dir"].as<std::string>()};

    // For now: Copy the options to the other tests
    ray_scan_cfg.track_generator() = hel_scan_cfg.track_generator();
    str_nav_cfg.propagation() = hel_nav_cfg.propagation();
    str_nav_cfg.fail_on_diff(false);

    detector_t::geometry_context ctx{};
    vecmem::host_memory_resource host_mr;

    const auto [det, names] =
        detray::io::read_detector<detector_t>(host_mr, reader_cfg);
    const std::string& det_name = det.name(names);

    // Create the whiteboard for data transfer between the steps
    auto white_board = std::make_shared<test::whiteboard>();
    const std::string file_prefix{data_dir + "/" + det_name};
    ray_scan_cfg.name(det_name + "_ray_scan");
    ray_scan_cfg.intersection_file(file_prefix + "_ray_scan_intersections");
    ray_scan_cfg.track_param_file(file_prefix + "_ray_scan_track_parameters");

    hel_scan_cfg.name(det_name + "_helix_scan");
    // Let the Newton algorithm dynamically choose tol. based on approx. error
    hel_scan_cfg.mask_tolerance({detray::detail::invalid_value<scalar>(),
                                 detray::detail::invalid_value<scalar>()});
    hel_scan_cfg.intersection_file(file_prefix + "_helix_scan_intersections");
    hel_scan_cfg.track_param_file(file_prefix + "_helix_scan_track_parameters");

    // General data consistency of the detector
    detray::test::register_checks<detray::test::consistency_check>(
        det, names, con_chk_cfg, ctx);

    // Navigation link consistency, discovered by ray intersection
    detray::test::register_checks<detray::test::ray_scan>(
        det, names, ray_scan_cfg, ctx, white_board);

    // Comparison of straight line navigation with ray scan
    str_nav_cfg.name(det_name + "_straight_line_navigation");
    // Number of tracks to check
    str_nav_cfg.n_tracks(ray_scan_cfg.track_generator().n_tracks());
    // Ensure that the same mask tolerance is used
    auto mask_tolerance = ray_scan_cfg.mask_tolerance();
    str_nav_cfg.propagation().navigation.min_mask_tolerance =
        static_cast<float>(mask_tolerance[0]);
    str_nav_cfg.propagation().navigation.max_mask_tolerance =
        static_cast<float>(mask_tolerance[1]);
    str_nav_cfg.intersection_file(ray_scan_cfg.intersection_file());
    str_nav_cfg.track_param_file(ray_scan_cfg.track_param_file());

    detray::test::register_checks<detray::test::straight_line_navigation>(
        det, names, str_nav_cfg, ctx, white_board);

    // Navigation link consistency, discovered by helix intersection
    detray::test::register_checks<detray::test::helix_scan>(
        det, names, hel_scan_cfg, ctx, white_board);

    // Comparison of navigation in a constant B-field with helix
    hel_nav_cfg.name(det_name + "_helix_navigation");
    hel_nav_cfg.fail_on_diff(false);
    // Number of tracks to check
    hel_nav_cfg.n_tracks(hel_scan_cfg.track_generator().n_tracks());
    hel_nav_cfg.p_range(hel_scan_cfg.track_generator().mom_range());
    hel_nav_cfg.intersection_file(hel_scan_cfg.intersection_file());
    hel_nav_cfg.track_param_file(hel_scan_cfg.track_param_file());

    detray::test::register_checks<detray::test::helix_navigation>(
        det, names, hel_nav_cfg, ctx, white_board);

    // Run the checks
    return RUN_ALL_TESTS();
}
