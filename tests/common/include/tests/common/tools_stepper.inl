/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/transform_store.hpp"
#include "detray/core/type_registry.hpp"
#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/masks/rectangle2.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/track_generators.hpp"

// google-test include(s)
#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;
using vector2 = __plugin::vector2<scalar>;
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;
using transform3 = __plugin::transform3<scalar>;
using matrix_operator = standard_matrix_operator<scalar>;
using mag_field_t = constant_magnetic_field<>;
using rk_stepper_t = rk_stepper<mag_field_t, transform3>;
using crk_stepper_t = rk_stepper<mag_field_t, transform3, constrained_step<>>;

namespace {

constexpr scalar epsilon = 1e-3;
constexpr scalar path_limit = 100 * unit_constants::cm;

// dummy navigation struct
struct nav_state {
    scalar operator()() const { return _step_size; }
    inline auto current_object() const -> dindex { return dindex_invalid; }

    inline void set_full_trust() {}
    inline void set_high_trust() {}
    inline void set_fair_trust() {}
    inline void set_no_trust() {}
    inline bool abort() { return false; }

    scalar _step_size = 1. * unit_constants::mm;
};

// dummy propagator state
template <typename stepping_t, typename navigation_t>
struct prop_state {
    stepping_t _stepping;
    navigation_t _navigation;
};

}  // namespace

// This tests the base functionality of the line stepper
TEST(ALGEBRA_PLUGIN, line_stepper) {
    using namespace step;

    // Line stepper with and without constrained stepping
    using line_stepper_t = line_stepper<transform3>;
    using cline_stepper_t = line_stepper<transform3, constrained_step<>>;

    point3 pos{0., 0., 0.};
    vector3 mom{1., 1., 0.};
    free_track_parameters<transform3> track(pos, 0, mom, -1);
    free_track_parameters<transform3> c_track(pos, 0, mom, -1);

    line_stepper_t l_stepper;
    cline_stepper_t cl_stepper;

    prop_state<line_stepper_t::state, nav_state> propagation{
        line_stepper_t::state{track}, nav_state{}};
    prop_state<cline_stepper_t::state, nav_state> c_propagation{
        cline_stepper_t::state{c_track}, nav_state{}};

    cline_stepper_t::state &cl_state = c_propagation._stepping;

    // Test the setting of step constraints
    cl_state.template set_constraint<constraint::e_accuracy>(
        10 * unit_constants::mm);
    cl_state.template set_constraint<constraint::e_actor>(2 *
                                                          unit_constants::mm);
    cl_state.template set_constraint<constraint::e_aborter>(5 *
                                                            unit_constants::mm);
    cl_state.template set_constraint<constraint::e_user>(0.5 *
                                                         unit_constants::mm);
    ASSERT_NEAR(cl_state.constraints().template size<>(),
                0.5 * unit_constants::mm, epsilon);

    // Release all except 'actor', then set 'user' again
    cl_state.template release_step<constraint::e_accuracy>();
    cl_state.template release_step<constraint::e_aborter>();
    cl_state.template release_step<constraint::e_user>();
    ASSERT_NEAR(cl_state.constraints().template size<>(),
                2 * unit_constants::mm, epsilon);

    cl_state.template set_constraint<constraint::e_user>(0.5 *
                                                         unit_constants::mm);
    ASSERT_NEAR(cl_state.constraints().template size<>(),
                0.5 * unit_constants::mm, epsilon);

    // Run a few steps
    ASSERT_TRUE(l_stepper.step(propagation));
    // Step constraint to half step size
    ASSERT_TRUE(cl_stepper.step(c_propagation));
    ASSERT_TRUE(cl_stepper.step(c_propagation));

    track = propagation._stepping();
    ASSERT_FLOAT_EQ(track.pos()[0], 1. / std::sqrt(2));
    ASSERT_FLOAT_EQ(track.pos()[1], 1. / std::sqrt(2));
    ASSERT_FLOAT_EQ(track.pos()[2], 0.);

    c_track = c_propagation._stepping();
    ASSERT_FLOAT_EQ(c_track.pos()[0], 1. / std::sqrt(2));
    ASSERT_FLOAT_EQ(c_track.pos()[1], 1. / std::sqrt(2));
    ASSERT_FLOAT_EQ(c_track.pos()[2], 0.);

    ASSERT_TRUE(l_stepper.step(propagation));

    track = propagation._stepping();
    ASSERT_FLOAT_EQ(track.pos()[0], std::sqrt(2));
    ASSERT_FLOAT_EQ(track.pos()[1], std::sqrt(2));
    ASSERT_FLOAT_EQ(track.pos()[2], 0.);
}

// This tests the base functionality of the Runge-Kutta stepper
TEST(ALGEBRA_PLUGIN, rk_stepper) {
    using namespace step;

    // RK stepper configurations
    constexpr unsigned int theta_steps = 100;
    constexpr unsigned int phi_steps = 100;
    constexpr unsigned int rk_steps = 100;

    // Constant magnetic field
    vector3 B{1 * unit_constants::T, 1 * unit_constants::T,
              1 * unit_constants::T};
    mag_field_t mag_field(B);

    // RK stepper
    rk_stepper_t rk_stepper;
    crk_stepper_t crk_stepper;

    // Set origin position of tracks
    const point3 ori{0., 0., 0.};
    const scalar p_mag = 10;

    // Iterate through uniformly distributed momentum directions
    for (auto track :
         uniform_track_generator<free_track_parameters<transform3>>(
             theta_steps, phi_steps, ori, p_mag)) {
        // Generate track state used for propagation with constrained step size
        free_track_parameters c_track(track);

        // helix gun
        detail::helix helix(track, &B);

        // RK Stepping into forward direction
        prop_state<rk_stepper_t::state, nav_state> propagation{
            rk_stepper_t::state{track, mag_field}, nav_state{}};
        prop_state<crk_stepper_t::state, nav_state> c_propagation{
            crk_stepper_t::state{c_track, mag_field}, nav_state{}};

        rk_stepper_t::state &rk_state = propagation._stepping;
        crk_stepper_t::state &crk_state = c_propagation._stepping;

        // Retrieve one of the navigation states
        nav_state &n_state = propagation._navigation;
        nav_state &cn_state = c_propagation._navigation;

        crk_state.template set_constraint<constraint::e_user>(
            0.5 * unit_constants::mm);
        n_state._step_size = 1. * unit_constants::mm;
        cn_state._step_size = 1. * unit_constants::mm;
        ASSERT_NEAR(crk_state.constraints().template size<>(),
                    0.5 * unit_constants::mm, epsilon);

        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk_stepper.step(propagation);
            crk_stepper.step(c_propagation);
            crk_stepper.step(c_propagation);
        }

        // get relative error by dividing error with path length
        ASSERT_NEAR(rk_state.path_length(), crk_state.path_length(), epsilon);
        ASSERT_NEAR(getter::norm(rk_state().pos() - crk_state().pos()) /
                        rk_state.path_length(),
                    0, epsilon);

        const auto helix_pos = helix(rk_state.path_length());
        const auto forward_pos = rk_state().pos();
        const point3 forward_relative_error{(1. / rk_state.path_length()) *
                                            (forward_pos - helix_pos)};

        // Make sure that relative error is smaller than epsion
        EXPECT_NEAR(getter::norm(forward_relative_error), 0, epsilon);

        // Roll the same track back to the origin
        // Use the same path length, since there is no overstepping
        const scalar path_length = rk_state.path_length();
        n_state._step_size *= -1. * unit_constants::mm;
        cn_state._step_size *= -1. * unit_constants::mm;
        for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
            rk_stepper.step(propagation);
            crk_stepper.step(c_propagation);
            crk_stepper.step(c_propagation);
        }

        ASSERT_NEAR(rk_state.path_length(), crk_state.path_length(), epsilon);
        ASSERT_NEAR(getter::norm(rk_state().pos() - crk_state().pos()) /
                        (2 * path_length),
                    0, epsilon);

        const auto backward_pos = rk_state().pos();
        const point3 backward_relative_error{1. / (2. * path_length) *
                                             (backward_pos - ori)};

        // Make sure that relative error is smaller than epsion
        EXPECT_NEAR(getter::norm(backward_relative_error), 0, epsilon);
    }
}