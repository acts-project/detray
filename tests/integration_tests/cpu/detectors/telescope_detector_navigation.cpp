/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/build_telescope_detector.hpp"
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

    //
    // Telescope detector configuration
    //
    using tel_detector_t = detector<telescope_metadata<rectangle2D>>;
    using scalar_t = typename tel_detector_t::scalar_type;

    tel_det_config<rectangle2D> tel_cfg{20.f * unit<scalar_t>::mm,
                                        20.f * unit<scalar_t>::mm};
    tel_cfg.n_surfaces(10u).length(500.f * unit<scalar_t>::mm);

    vecmem::host_memory_resource host_mr;

    const auto [tel_det, tel_names] =
        build_telescope_detector(host_mr, tel_cfg);

    auto white_board = std::make_shared<test::whiteboard>();

    // General data consistency of the detector
    test::consistency_check<tel_detector_t>::config cfg_cons{};
    detail::register_checks<test::consistency_check>(
        tel_det, tel_names, cfg_cons.name("telescope_detector_consistency"));

    // Navigation link consistency, discovered by ray intersection
    test::ray_scan<tel_detector_t>::config cfg_ray_scan{};
    cfg_ray_scan.name("telescope_detector_ray_scan");
    cfg_ray_scan.whiteboard(white_board);
    cfg_ray_scan.track_generator().theta_steps(100u).phi_steps(100u);
    cfg_ray_scan.track_generator().origin({0.f, 0.f, -0.05f});
    cfg_ray_scan.track_generator().theta_range(constant<scalar_t>::pi_4,
                                               constant<scalar_t>::pi_2);

    detail::register_checks<test::ray_scan>(tel_det, tel_names, cfg_ray_scan);

    // Navigation link consistency, discovered by helix intersection
    test::helix_scan<tel_detector_t>::config cfg_hel_scan{};
    cfg_hel_scan.name("telescope_detector_helix_scan");
    cfg_hel_scan.whiteboard(white_board);
    cfg_hel_scan.track_generator().p_tot(10.f * unit<scalar_t>::GeV);
    cfg_hel_scan.track_generator().origin({0.f, 0.f, -0.05f});
    cfg_hel_scan.track_generator().theta_steps(100u).phi_steps(100u);
    cfg_hel_scan.track_generator().theta_range(constant<scalar_t>::pi_4,
                                               constant<scalar_t>::pi_2);
    detail::register_checks<test::helix_scan>(tel_det, tel_names, cfg_hel_scan);

    // Comparision of straight line navigation with ray scan
    test::straight_line_navigation<tel_detector_t>::config cfg_str_nav{};
    cfg_str_nav.name("telescope_detector_straight_line_navigation");
    cfg_str_nav.whiteboard(white_board);

    detail::register_checks<test::straight_line_navigation>(tel_det, tel_names,
                                                            cfg_str_nav);

    // Comparision of navigation in a constant B-field with helix
    test::helix_navigation<tel_detector_t>::config cfg_hel_nav{};
    cfg_hel_nav.name("telescope_detector_helix_navigation");
    cfg_hel_nav.whiteboard(white_board);

    detail::register_checks<test::helix_navigation>(tel_det, tel_names,
                                                    cfg_hel_nav);

    // Run the checks
    return RUN_ALL_TESTS();
}
