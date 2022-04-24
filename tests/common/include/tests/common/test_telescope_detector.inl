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

using vector3 = __plugin::vector3<detray::scalar>;
}

// This tests the construction and general methods of the navigator
TEST(ALGEBRA_PLUGIN, telescope_detector) {
    using namespace detray;

    using b_field_t = constant_magnetic_field<scalar>;
    using ln_stepper_t = line_stepper<free_track_parameters>;
    using rk_stepper_t = rk_stepper<b_field_t, free_track_parameters>;
    using inspector_t = navigation::print_inspector;

    // Test tolerance
    constexpr scalar tol = 1e-4;

    vecmem::host_memory_resource host_mr;

    vector3 B_z{0., 0., 1. * unit_constants::T};
    vector3 B_x{1. * unit_constants::T, 0., 0.};
    b_field_t b_field_z{B_z};
    b_field_t b_field_x{B_x};

    rk_stepper_t rk_step_z{b_field_z};
    rk_stepper_t rk_step_x{b_field_x};
    ln_stepper_t ln_step{};

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

    // Same telescope, but in x direction
    point3 pos{0., 0., 0.};
    vector3 mom{1., 0., 0.};
    free_track_parameters pilot_track(pos, 0, mom, -1);

    const auto x_tel_det = create_telescope_detector<rectangular>(
        host_mr, n_surfaces, tel_length, pilot_track, ln_step);

    // Telescope navigation should be symmetric in x and z
    pos = {0., 0., 0.};
    mom = {0., 0., 1.};
    free_track_parameters test_track1(pos, 0, mom, -1);
    free_track_parameters test_track2(pos, 0, mom, -1);
    mom = {1., 0., 0.};
    free_track_parameters test_track3(pos, 0, mom, -1);

    rk_stepper_t::state s_z1(test_track1);
    rk_stepper_t::state s_z2(test_track2);
    rk_stepper_t::state s_x(test_track3);

    navigator<decltype(z_tel_det1), inspector_t> nav_z1(z_tel_det1);
    navigator<decltype(z_tel_det2), inspector_t> nav_z2(z_tel_det2);
    navigator<decltype(x_tel_det), inspector_t> nav_x(x_tel_det);
    decltype(nav_z1)::state n_z1, n_z2, n_x;

    bool heartbeat_z1 = nav_z1.init(n_z1, s_z1);
    bool heartbeat_z2 = nav_z2.init(n_z2, s_z2);
    bool heartbeat_x = nav_x.init(n_x, s_x);

    while (heartbeat_z1 and heartbeat_z2 and heartbeat_x) {

        EXPECT_TRUE(heartbeat_z1);
        EXPECT_TRUE(heartbeat_z2);
        EXPECT_TRUE(heartbeat_x);

        heartbeat_z1 &= rk_step_x.step(s_z1, n_z1);
        heartbeat_z2 &= rk_step_x.step(s_z2, n_z2);
        heartbeat_x &= rk_step_z.step(s_x, n_x);

        heartbeat_z1 &= nav_z1.update(n_z1, s_z1);
        heartbeat_z2 &= nav_z2.update(n_z2, s_z2);
        heartbeat_x &= nav_x.update(n_x, s_x);

        EXPECT_NEAR(std::fabs(s_z1._path_length - s_z2._path_length) /
                        s_z1._path_length,
                    0., tol);
        EXPECT_NEAR(
            std::fabs(s_z1._path_length - s_x._path_length) / s_x._path_length,
            0., tol);
        EXPECT_NEAR(getter::norm(test_track1.pos() - test_track2.pos()) /
                        getter::norm(test_track1.pos()),
                    0., tol);
        EXPECT_NEAR(getter::norm(test_track1.dir() - test_track2.dir()) /
                        getter::norm(test_track1.dir()),
                    0., tol);
    }
    ASSERT_TRUE(n_z1.is_complete()) << n_z1.inspector().to_string();
    ASSERT_TRUE(n_z2.is_complete()) << n_z2.inspector().to_string();
    ASSERT_TRUE(n_x.is_complete()) << n_x.inspector().to_string();

    //
    // Build a telescope along a bent track
    //
    pos = {0., 0., 0.};
    mom = {0., 1., 0.};
    pilot_track = free_track_parameters(pos, 0, mom, -1);
    pilot_track.set_overstep_tolerance(-10 * unit_constants::um);

    const auto tel_det = create_telescope_detector<rectangular>(
        host_mr, n_surfaces, tel_length, pilot_track, rk_step_z);

    // make at least sure it is navigatable
    rk_stepper_t::state s_tel(pilot_track);
    navigator<decltype(tel_det), inspector_t> nav_tel(tel_det);
    decltype(nav_tel)::state n_tel;
    bool heartbeat_tel = nav_tel.init(n_tel, s_tel);

    while (heartbeat_tel) {
        heartbeat_tel &= rk_step_z.step(s_tel, n_tel);
        heartbeat_tel &= nav_tel.update(n_tel, s_tel);
    }
    ASSERT_TRUE(n_tel.is_complete()) << n_tel.inspector().to_string();
}
