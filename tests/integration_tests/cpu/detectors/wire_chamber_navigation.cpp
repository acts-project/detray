/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_wire_chamber.hpp"
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
    // Wire Chamber configuration
    //
    vecmem::host_memory_resource host_mr;

    using wire_chamber_t = detector<>;
    using scalar_t = typename wire_chamber_t::scalar_type;

    wire_chamber_config wire_chamber_cfg{};
    wire_chamber_cfg.half_z(500.f * unit<scalar>::mm);

    auto [det, names] = create_wire_chamber(host_mr, wire_chamber_cfg);

    auto white_board = std::make_shared<test::whiteboard>();

    // General data consistency of the detector
    test::consistency_check<wire_chamber_t>::config cfg_cons{};
    detail::register_checks<test::consistency_check>(
        det, names, cfg_cons.name("wire_chamber_consistency"));

    // Navigation link consistency, discovered by ray intersection
    test::ray_scan<wire_chamber_t>::config cfg_ray_scan{};
    cfg_ray_scan.name("wire_chamber_ray_scan");
    cfg_ray_scan.whiteboard(white_board);
    cfg_ray_scan.track_generator().n_tracks(10000u);

    detail::register_checks<test::ray_scan>(det, names, cfg_ray_scan);

    // Navigation link consistency, discovered by helix intersection
    test::helix_scan<wire_chamber_t>::config cfg_hel_scan{};
    cfg_hel_scan.name("wire_chamber_helix_scan");
    cfg_hel_scan.whiteboard(white_board);
    cfg_hel_scan.track_generator().n_tracks(10000u);
    // TODO: Fails for smaller momenta
    cfg_hel_scan.track_generator().p_T(3.f * unit<scalar_t>::GeV);

    detail::register_checks<test::helix_scan>(det, names, cfg_hel_scan);

    // Comparison of straight line navigation with ray scan
    test::straight_line_navigation<wire_chamber_t>::config cfg_str_nav{};
    cfg_str_nav.name("wire_chamber_straight_line_navigation");
    cfg_str_nav.whiteboard(white_board);
    cfg_str_nav.propagation().navigation.search_window = {2u, 2u};
    auto mask_tolerance = cfg_ray_scan.mask_tolerance();
    cfg_str_nav.propagation().navigation.min_mask_tolerance = mask_tolerance[0];
    cfg_str_nav.propagation().navigation.max_mask_tolerance = mask_tolerance[1];

    detail::register_checks<test::straight_line_navigation>(det, names,
                                                            cfg_str_nav);

    // Comparison of navigation in a constant B-field with helix
    test::helix_navigation<wire_chamber_t>::config cfg_hel_nav{};
    cfg_hel_nav.name("wire_chamber_helix_navigation");
    cfg_hel_nav.whiteboard(white_board);
    cfg_hel_nav.propagation().navigation.search_window = {3u, 3u};

    detail::register_checks<test::helix_navigation>(det, names, cfg_hel_nav);

    // Run the checks
    return RUN_ALL_TESTS();
}
