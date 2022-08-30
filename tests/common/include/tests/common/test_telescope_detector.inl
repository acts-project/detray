/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/propagator/track.hpp"
#include "tests/common/tools/create_telescope_detector.hpp"
#include "tests/common/tools/inspectors.hpp"

/// @note __plugin has to be defined with a preprocessor command
namespace detray {

namespace {

using vector3 = __plugin::vector3<detray::scalar>;

// dummy propagator state
template <typename stepping_t, typename navigation_t>
struct prop_state {

    stepping_t _stepping;
    navigation_t _navigation;
    using field_type = typename stepping_t::field_type;

    template <typename track_t>
    prop_state(const track_t &t_in, const field_type &field,
               const typename navigation_t::detector_type &det)
        : _stepping(t_in, field), _navigation(det) {}
};

}  // anonymous namespace

}  // namespace detray

// This tests the construction and general methods of the navigator
TEST(ALGEBRA_PLUGIN, telescope_detector) {

    using namespace detray;

    using b_field_t = constant_magnetic_field<>;
    using ln_stepper_t = line_stepper<free_track_parameters>;
    using rk_stepper_t = rk_stepper<b_field_t, free_track_parameters>;
    using inspector_t = navigation::print_inspector;

    // Use rectangular surfaces
    constexpr bool rectangular = false;

    // Test tolerance
    constexpr scalar tol = 1e-4;

    vecmem::host_memory_resource host_mr;

    // B-fields
    vector3 B_z{0., 0., 1. * unit_constants::T};
    vector3 B_x{1. * unit_constants::T, 0., 0.};
    b_field_t b_field_z{B_z};
    b_field_t b_field_x{B_x};

    // steppers
    rk_stepper_t rk_stepper_z;
    rk_stepper_t rk_stepper_x;
    ln_stepper_t ln_stepper;

    //
    // telescope along z
    //

    // Build from given module positions
    std::vector<scalar> positions = {0.,   50., 100., 150., 200., 250.,
                                     300., 350, 400,  450., 500.};
    // Build telescope detector with unbounded planes
    const auto z_tel_det1 =
        create_telescope_detector<rectangular>(host_mr, positions);

    // Build the same telescope detector with rectangular planes and given
    // length/number of surfaces
    dindex n_surfaces = 11;
    scalar tel_length = 500. * unit_constants::mm;
    const auto z_tel_det2 =
        create_telescope_detector<rectangular>(host_mr, n_surfaces, tel_length);

    // Compare
    for (std::size_t i = 0; i < z_tel_det1.surfaces().size(); ++i) {
        EXPECT_TRUE(z_tel_det1.surface_by_index(i) ==
                    z_tel_det2.surface_by_index(i));
    }

    //
    // telescope along x
    //

    // Same telescope, but in x direction and created from custom stepper
    point3 pos{0., 0., 0.};
    vector3 mom{1., 0., 0.};
    free_track_parameters pilot_track(pos, 0, mom, -1);
    typename ln_stepper_t::state ln_stepping(pilot_track);

    const auto x_tel_det = create_telescope_detector<rectangular>(
        host_mr, n_surfaces, tel_length, ln_stepper, ln_stepping);

    //
    // test propagation in all telescope detector instances
    //

    // Telescope navigation should be symmetric in x and z
    pos = {0., 0., 0.};
    mom = {0., 0., 1.};
    free_track_parameters test_track_z1(pos, 0, mom, -1);
    free_track_parameters test_track_z2(pos, 0, mom, -1);
    mom = {1., 0., 0.};
    free_track_parameters test_track_x(pos, 0, mom, -1);

    // navigators
    navigator<decltype(z_tel_det1), inspector_t> navigator_z1;
    navigator<decltype(z_tel_det2), inspector_t> navigator_z2;
    navigator<decltype(x_tel_det), inspector_t> navigator_x;
    using navigation_state_t = decltype(navigator_z1)::state;
    using stepping_state_t = rk_stepper_t::state;

    // propagation states
    prop_state<stepping_state_t, navigation_state_t> propgation_z1(
        test_track_z1, b_field_z, z_tel_det1);
    prop_state<stepping_state_t, navigation_state_t> propgation_z2(
        test_track_z2, b_field_z, z_tel_det2);
    prop_state<stepping_state_t, navigation_state_t> propgation_x(
        test_track_x, b_field_x, x_tel_det);

    stepping_state_t &stepping_z1 = propgation_z1._stepping;
    stepping_state_t &stepping_z2 = propgation_z2._stepping;
    stepping_state_t &stepping_x = propgation_x._stepping;

    navigation_state_t &navigation_z1 = propgation_z1._navigation;
    navigation_state_t &navigation_z2 = propgation_z2._navigation;
    navigation_state_t &navigation_x = propgation_x._navigation;

    // propagate all telescopes
    bool heartbeat_z1 = navigator_z1.init(propgation_z1);
    bool heartbeat_z2 = navigator_z2.init(propgation_z2);
    bool heartbeat_x = navigator_x.init(propgation_x);

    while (heartbeat_z1 and heartbeat_z2 and heartbeat_x) {

        // check that all propagation flows are still running
        EXPECT_TRUE(heartbeat_z1);
        EXPECT_TRUE(heartbeat_z2);
        EXPECT_TRUE(heartbeat_x);

        heartbeat_z1 &= rk_stepper_z.step(propgation_z1);
        heartbeat_z2 &= rk_stepper_z.step(propgation_z2);
        heartbeat_x &= rk_stepper_x.step(propgation_x);

        navigation_z1.set_high_trust();
        navigation_z2.set_high_trust();
        navigation_x.set_high_trust();

        heartbeat_z1 &= navigator_z1.update(propgation_z1);
        heartbeat_z2 &= navigator_z2.update(propgation_z2);
        heartbeat_x &= navigator_x.update(propgation_x);
        // The track path lengths should match between all propagations
        EXPECT_NEAR(
            std::fabs(stepping_z1._path_length - stepping_z2._path_length) /
                stepping_z1._path_length,
            0., tol);
        EXPECT_NEAR(
            std::fabs(stepping_z1._path_length - stepping_x._path_length) /
                stepping_x._path_length,
            0., tol);
        // The track positions in z should match exactly
        EXPECT_NEAR(getter::norm(stepping_z1().pos() - stepping_z2().pos()) /
                        getter::norm(stepping_z1().pos()),
                    0., tol);
        EXPECT_NEAR(getter::norm(stepping_z1().dir() - stepping_z2().dir()) /
                        getter::norm(stepping_z1().dir()),
                    0., tol);
    }

    // check that all propagation flows exited successfully
    ASSERT_TRUE(navigation_z1.is_complete())
        << navigation_z1.inspector().to_string();
    ASSERT_TRUE(navigation_z2.is_complete())
        << navigation_z2.inspector().to_string();
    ASSERT_TRUE(navigation_x.is_complete())
        << navigation_x.inspector().to_string();

    //
    // Build a telescope along a bent track
    //
    pos = {0., 0., 0.};
    mom = {0., 1., 0.};
    pilot_track = free_track_parameters(pos, 0, mom, -1);
    pilot_track.set_overstep_tolerance(-10 * unit_constants::um);

    typename rk_stepper_t::state rk_stepping_z(pilot_track, b_field_z);
    typename rk_stepper_t::state rk_stepping_x(pilot_track, b_field_x);

    const auto tel_detector = create_telescope_detector<rectangular>(
        host_mr, n_surfaces, tel_length, rk_stepper_z, rk_stepping_z);

    // make at least sure it is navigatable
    navigator<decltype(tel_detector), inspector_t> tel_navigator;

    prop_state<stepping_state_t, navigation_state_t> tel_propagation(
        pilot_track, b_field_z, tel_detector);
    navigation_state_t &tel_navigation = tel_propagation._navigation;

    // run propagation
    bool heartbeat_tel = tel_navigator.init(tel_propagation);

    while (heartbeat_tel) {
        heartbeat_tel &= rk_stepper_z.step(tel_propagation);
        tel_navigation.set_high_trust();
        heartbeat_tel &= tel_navigator.update(tel_propagation);
    }
    // check that propagation was successful
    ASSERT_TRUE(tel_navigation.is_complete())
        << tel_navigation.inspector().to_string();
}
