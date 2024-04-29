/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/build_toy_detector.hpp"
#include "detray/test/detail/register_checks.hpp"
#include "detray/test/detail/whiteboard.hpp"
#include "detray/test/detector_consistency.hpp"
#include "detray/test/detector_helix_scan.hpp"
#include "detray/test/detector_ray_scan.hpp"
#include "detray/test/helix_navigation.hpp"
#include "detray/test/straight_line_navigation.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

int main(int argc, char **argv) {

    using namespace detray;

    // Filter out the google test flags
    ::testing::InitGoogleTest(&argc, argv);

    using toy_detector_t = detector<toy_metadata>;
    using scalar_t = typename toy_detector_t::scalar_type;

    //
    // Toy detector configuration
    //
    toy_det_config<scalar_t> toy_cfg{};
    toy_cfg.n_brl_layers(4u).n_edc_layers(7u);

    // Build the geometry
    vecmem::host_memory_resource host_mr;
    auto [toy_det, toy_names] = build_toy_detector(host_mr, toy_cfg);

    auto white_board = std::make_shared<test::whiteboard>();

    // General data consistency of the detector
    test::consistency_check<toy_detector_t>::config cfg_cons{};
    detail::register_checks<test::consistency_check>(
        toy_det, toy_names, cfg_cons.name("toy_detector_consistency"));

    // Navigation link consistency, discovered by ray intersection
    test::ray_scan<toy_detector_t>::config cfg_ray_scan{};
    cfg_ray_scan.name("toy_detector_ray_scan");
    cfg_ray_scan.whiteboard(white_board);
    cfg_ray_scan.track_generator().n_tracks(10000u);

    detail::register_checks<test::ray_scan>(toy_det, toy_names, cfg_ray_scan);

    // Navigation link consistency, discovered by helix intersection
    test::helix_scan<toy_detector_t>::config cfg_hel_scan{};
    cfg_hel_scan.name("toy_detector_helix_scan");
    cfg_hel_scan.whiteboard(white_board);
    cfg_hel_scan.track_generator().n_tracks(10000u);
    cfg_hel_scan.track_generator().eta_range(-4.f, 4.f);
    cfg_hel_scan.track_generator().p_T(1.f * unit<scalar_t>::GeV);

    detail::register_checks<test::helix_scan>(toy_det, toy_names, cfg_hel_scan);

    // Comparison of straight line navigation with ray scan
    test::straight_line_navigation<toy_detector_t>::config cfg_str_nav{};
    cfg_str_nav.name("toy_detector_straight_line_navigation");
    cfg_str_nav.whiteboard(white_board);
    cfg_str_nav.propagation().navigation.search_window = {3u, 3u};
    auto mask_tolerance = cfg_ray_scan.mask_tolerance();
    cfg_str_nav.propagation().navigation.min_mask_tolerance = mask_tolerance[0];
    cfg_str_nav.propagation().navigation.max_mask_tolerance = mask_tolerance[1];

    detail::register_checks<test::straight_line_navigation>(toy_det, toy_names,
                                                            cfg_str_nav);

    // Comparison of navigation in a constant B-field with helix
    test::helix_navigation<toy_detector_t>::config cfg_hel_nav{};
    cfg_hel_nav.name("toy_detector_helix_navigation");
    cfg_hel_nav.whiteboard(white_board);
    cfg_hel_nav.propagation().navigation.search_window = {3u, 3u};

    detail::register_checks<test::helix_navigation>(toy_det, toy_names,
                                                    cfg_hel_nav);

    // Run the checks
    return RUN_ALL_TESTS();
}
