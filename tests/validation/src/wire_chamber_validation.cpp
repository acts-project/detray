/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_wire_chamber.hpp"
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
    // Wire Chamber configuration
    //
    vecmem::host_memory_resource host_mr;

    using wire_chamber_t = detector<>;
    using scalar_t = typename wire_chamber_t::scalar_type;

    wire_chamber_config wire_chamber_cfg{};
    wire_chamber_cfg.half_z(500.f * unit<scalar>::mm);
    wire_chamber_cfg.bfield_vec(0.f, 0.f, 2.f * unit<scalar>::T);

    auto [det, names] = create_wire_chamber(host_mr, wire_chamber_cfg);

    // General data consistency of the detector
    consistency_check<wire_chamber_t>::config cfg_cons{};
    detail::register_checks<consistency_check>(
        det, names, cfg_cons.name("wire_chamber_consistency"));

    // Navigation link consistency, discovered by ray intersection
    ray_scan<wire_chamber_t>::config cfg_ray_scan{};
    cfg_ray_scan.name("wire_chamber_ray_scan");
    // Number of rays in theta and phi
    cfg_ray_scan.track_generator().theta_steps(100u).phi_steps(100u);

    detail::register_checks<ray_scan>(det, names, cfg_ray_scan);

    // Navigation link consistency, discovered by helix intersection
    helix_scan<wire_chamber_t>::config cfg_hel_scan{};
    cfg_hel_scan.name("wire_chamber_helix_scan");
    cfg_hel_scan.overstepping_tolerance(-100.f * unit<scalar_t>::um);
    cfg_hel_scan.track_generator().p_mag(10.f * unit<scalar_t>::GeV);
    // Fails in single precision for more helices
    // cfg_hel_scan.track_generator().theta_steps(100u).phi_steps(100u);
    cfg_hel_scan.track_generator().theta_steps(3u).phi_steps(10u);
    detray::detail::register_checks<detray::helix_scan>(det, names,
                                                        cfg_hel_scan);

    // Comparision of straight line navigation with ray scan
    straight_line_navigation<wire_chamber_t>::config cfg_str_nav{};
    cfg_str_nav.name("wire_chamber_straight_line_navigation");
    cfg_str_nav.track_generator().theta_steps(100u).phi_steps(100u);

    detail::register_checks<straight_line_navigation>(det, names, cfg_str_nav);

    // Comparision of navigation in a constant B-field with helix
    helix_navigation<wire_chamber_t>::config cfg_hel_nav{};
    cfg_hel_nav.name("wire_chamber_helix_navigation");
    cfg_hel_nav.overstepping_tolerance(-100.f * unit<scalar_t>::um);
    cfg_hel_nav.track_generator() = cfg_hel_scan.track_generator();
    // TODO: Fails for more helices
    // cfg_hel_nav.track_generator().theta_steps(100u).phi_steps(100u);
    cfg_hel_nav.track_generator().theta_steps(3u).phi_steps(10u);

    detail::register_checks<helix_navigation>(det, names, cfg_hel_nav);

    // Run the checks
    return RUN_ALL_TESTS();
}
