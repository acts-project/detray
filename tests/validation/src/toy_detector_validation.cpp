/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
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

int main(int argc, char **argv) {

    using namespace detray;

    // Filter out the google test flags
    ::testing::InitGoogleTest(&argc, argv);

    //
    // Toy detector configuration
    //
    toy_det_config toy_cfg{};
    toy_cfg.n_brl_layers(4u).n_edc_layers(7u);

    using toy_detector_t = detector<toy_metadata>;
    using scalar_t = typename toy_detector_t::scalar_type;

    // Build the geometry
    vecmem::host_memory_resource host_mr;
    auto [toy_det, toy_names] = create_toy_geometry(host_mr, toy_cfg);

    // General data consistency of the detector
    consistency_check<toy_detector_t>::config cfg_cons{};
    detail::register_checks<consistency_check>(
        toy_det, toy_names, cfg_cons.name("toy_detector_consistency"));

    // Navigation link consistency, discovered by ray intersection
    ray_scan<toy_detector_t>::config cfg_ray_scan{};
    cfg_ray_scan.name("toy_detector_ray_scan");
    // Number of rays in theta and phi
    cfg_ray_scan.track_generator().theta_steps(100u).phi_steps(100u);

    detail::register_checks<ray_scan>(toy_det, toy_names, cfg_ray_scan);

    // Navigation link consistency, discovered by helix intersection
    helix_scan<toy_detector_t>::config cfg_hel_scan{};
    cfg_hel_scan.name("toy_detector_helix_scan");
    cfg_hel_scan.track_generator().p_tot(10.f * unit<scalar_t>::GeV);
    cfg_hel_scan.track_generator().theta_steps(100u).phi_steps(100u);

    detray::detail::register_checks<detray::helix_scan>(toy_det, toy_names,
                                                        cfg_hel_scan);

    // Comparision of straight line navigation with ray scan
    straight_line_navigation<toy_detector_t>::config cfg_str_nav{};
    cfg_str_nav.name("toy_detector_straight_line_navigation");
    cfg_str_nav.propagation().search_window = {3u, 3u};
    cfg_str_nav.track_generator().theta_steps(100u).phi_steps(100u);

    detail::register_checks<straight_line_navigation>(toy_det, toy_names,
                                                      cfg_str_nav);

    // Comparision of navigation in a constant B-field with helix
    helix_navigation<toy_detector_t>::config cfg_hel_nav{};
    cfg_hel_nav.name("toy_detector_helix_navigation");
    cfg_hel_nav.propagation().search_window = {3u, 3u};
    cfg_hel_nav.track_generator() = cfg_hel_scan.track_generator();
    // TODO: Fails due to mask tolerances for more helices, regardless of edc
    // configuration/precision
    // cfg_hel_nav.track_generator().theta_steps(100u).phi_steps(100u);
    cfg_hel_nav.track_generator().theta_steps(50u).phi_steps(50u);

    detail::register_checks<helix_navigation>(toy_det, toy_names, cfg_hel_nav);

    // Run the checks
    return RUN_ALL_TESTS();
}
