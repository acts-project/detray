/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/definitions/units.hpp"
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/navigation_policies.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/inspectors.hpp"

// This tests the construction and general methods of the navigator
TEST(ALGEBRA_PLUGIN, guided_navigator) {
    using namespace detray;
    using namespace navigation;
    using transform3_type = __plugin::transform3<scalar>;

    vecmem::host_memory_resource host_mr;

    // Use unbounded surfaces
    constexpr bool unbounded = true;

    detail::ray<transform3_type> default_trk({0, 0, 0}, 0, {0, 0, 1}, -1);

    // Module positions along z-axis
    const std::vector<scalar> positions = {0.,  10., 20., 30., 40., 50.,
                                           60., 70,  80,  90., 100.};
    // Build telescope detector with unbounded planes
    const auto telescope_det =
        create_telescope_detector<unbounded>(host_mr, positions, default_trk);

    // Inspectors are optional, of course
    using object_tracer_t =
        object_tracer<dvector, status::e_on_portal, status::e_on_module>;
    using inspector_t = aggregate_inspector<object_tracer_t, print_inspector>;
    using b_field_t = decltype(telescope_det)::bfield_type;
    using runge_kutta_stepper =
        rk_stepper<b_field_t::view_t, transform3_type, unconstrained_step,
                   guided_navigation>;
    using guided_navigator = navigator<decltype(telescope_det), inspector_t>;
    using actor_chain_t = actor_chain<dtuple, pathlimit_aborter>;
    using propagator_t =
        propagator<runge_kutta_stepper, guided_navigator, actor_chain_t>;

    // track must point into the direction of the telescope
    const point3 pos{0., 0., 0.};
    const vector3 mom{0., 0., 1.};
    free_track_parameters<transform3_type> track(pos, 0, mom, -1);
    const vector3 B{0, 0, 1 * unit<scalar>::T};
    const b_field_t b_field(
        b_field_t::backend_t::configuration_t{B[0], B[1], B[2]});

    // Actors
    pathlimit_aborter::state pathlimit{200. * unit<scalar>::cm};

    // Propagator
    propagator_t p(runge_kutta_stepper{}, guided_navigator{});
    propagator_t::state guided_state(track, b_field, telescope_det);

    // Propagate
    p.propagate(guided_state, std::tie(pathlimit));

    auto &nav_state = guided_state._navigation;
    auto &debug_printer = nav_state.inspector().template get<print_inspector>();
    auto &obj_tracer = nav_state.inspector().template get<object_tracer_t>();

    // Check that navigator exited
    ASSERT_TRUE(nav_state.is_complete()) << debug_printer.to_string();

    // Sequence of surface ids we expect to see
    const std::vector<dindex> sf_sequence = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    // Check the surfaces that have been visited by the navigation
    EXPECT_EQ(obj_tracer.object_trace.size(), sf_sequence.size());
    for (size_t i = 0; i < sf_sequence.size(); ++i) {
        const auto candidate = obj_tracer.object_trace[i];
        auto bcd = geometry::barcode{};
        bcd.set_volume(0UL).set_index(sf_sequence[i]);
        bcd.set_id((i == 0 or i == 10) ? surface_id::e_portal
                                       : surface_id::e_sensitive);
        EXPECT_TRUE(candidate.barcode == bcd)
            << "error at intersection on surface: " << candidate.barcode;
    }
}
