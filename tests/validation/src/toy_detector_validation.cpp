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
    constexpr std::size_t n_brl_layers{4u};
    constexpr std::size_t n_edc_layers{7u};
    vecmem::host_memory_resource host_mr;

    using toy_detector_t = detector<toy_metadata<>>;
    using b_field_t = toy_detector_t::bfield_type;
    using scalar_t = typename toy_detector_t::scalar_type;

    const toy_detector_t::vector3 B{
        0.f * unit<scalar>::T, 0.f * unit<scalar>::T, 2.f * unit<scalar>::T};

    const auto [toy_det, toy_names] = create_toy_geometry(
        host_mr,
        b_field_t(b_field_t::backend_t::configuration_t{B[0], B[1], B[2]}),
        n_brl_layers, n_edc_layers);

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
    cfg_hel_scan.overstepping_tolerance(-100.f * unit<scalar_t>::um);
    cfg_hel_scan.track_generator().p_mag(10.f * unit<scalar_t>::GeV);
    // Fails in single precision for more helices (unless is configured with
    // only three edc layers)
    // cfg_hel_scan.track_generator().theta_steps(100u).phi_steps(100u);
    cfg_hel_scan.track_generator().theta_steps(30u).phi_steps(30u);
    detray::detail::register_checks<detray::helix_scan>(toy_det, toy_names,
                                                        cfg_hel_scan);

    // Comparision of straight line navigation with ray scan
    straight_line_navigation<toy_detector_t>::config cfg_str_nav{};
    cfg_str_nav.name("toy_detector_straight_line_navigation");
    cfg_str_nav.track_generator().theta_steps(100u).phi_steps(100u);

    detail::register_checks<straight_line_navigation>(toy_det, toy_names,
                                                      cfg_str_nav);

    // Comparision of navigation in a constant B-field with helix
    helix_navigation<toy_detector_t>::config cfg_hel_nav{};
    cfg_hel_nav.name("toy_detector_helix_navigation");
    cfg_hel_nav.overstepping_tolerance(-100.f * unit<scalar_t>::um);
    cfg_hel_nav.track_generator() = cfg_hel_scan.track_generator();
    // TODO: Fails due to mask tolerances for more helices, regardless of edc
    // configuration/precision
    // cfg_hel_nav.track_generator().theta_steps(100u).phi_steps(100u);
    cfg_hel_nav.track_generator().theta_steps(30u).phi_steps(30u);

    detail::register_checks<helix_navigation>(toy_det, toy_names, cfg_hel_nav);

    // Run the checks
    return RUN_ALL_TESTS();
}
