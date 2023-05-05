/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// detray include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/tracks.hpp"

// google-test include(s)
#include <gtest/gtest.h>

// covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>
#include <covfie/core/vector.hpp>

using namespace detray;
using vector2 = test::vector2;
using vector3 = test::vector3;
using point3 = test::point3;
using transform3 = test::transform3;
using matrix_operator = standard_matrix_operator<scalar>;
using mag_field_t = covfie::field<covfie::backend::constant<
    covfie::vector::vector_d<scalar, 3>, covfie::vector::vector_d<scalar, 3>>>;
using rk_stepper_t = rk_stepper<mag_field_t::view_t, transform3>;
using crk_stepper_t =
    rk_stepper<mag_field_t::view_t, transform3, constrained_step<>>;

namespace {

constexpr scalar tol{1e-3f};

// dummy navigation struct
struct nav_state {
    scalar operator()() const { return _step_size; }
    inline auto current_object() const -> dindex { return dindex_invalid; }

    inline void set_full_trust() {}
    inline void set_high_trust() {}
    inline void set_fair_trust() {}
    inline void set_no_trust() {}
    inline bool abort() { return false; }

    scalar _step_size{1.f * unit<scalar>::mm};
};

// dummy propagator state
template <typename stepping_t, typename navigation_t>
struct prop_state {
    stepping_t _stepping;
    navigation_t _navigation;
};

}  // namespace

// This tests the base functionality of the line stepper
GTEST_TEST(detray_propagator, line_stepper) {
    using namespace step;

    // Line stepper with and without constrained stepping
    using line_stepper_t = line_stepper<transform3>;
    using cline_stepper_t = line_stepper<transform3, constrained_step<>>;

    point3 pos{0.f, 0.f, 0.f};
    vector3 mom{1.f, 1.f, 0.f};
    free_track_parameters<transform3> track(pos, 0.f, mom, -1.f);
    free_track_parameters<transform3> c_track(pos, 0.f, mom, -1.f);

    line_stepper_t l_stepper;
    cline_stepper_t cl_stepper;

    prop_state<line_stepper_t::state, nav_state> propagation{
        line_stepper_t::state{track}, nav_state{}};
    prop_state<cline_stepper_t::state, nav_state> c_propagation{
        cline_stepper_t::state{c_track}, nav_state{}};

    cline_stepper_t::state &cl_state = c_propagation._stepping;

    // Test the setting of step constraints
    cl_state.template set_constraint<constraint::e_accuracy>(10.f *
                                                             unit<scalar>::mm);
    cl_state.template set_constraint<constraint::e_actor>(2.f *
                                                          unit<scalar>::mm);
    cl_state.template set_constraint<constraint::e_aborter>(5.f *
                                                            unit<scalar>::mm);
    cl_state.template set_constraint<constraint::e_user>(0.5f *
                                                         unit<scalar>::mm);
    ASSERT_NEAR(cl_state.constraints().template size<>(),
                0.5f * unit<scalar>::mm, tol);

    // Release all except 'actor', then set 'user' again
    cl_state.template release_step<constraint::e_accuracy>();
    cl_state.template release_step<constraint::e_aborter>();
    cl_state.template release_step<constraint::e_user>();
    ASSERT_NEAR(cl_state.constraints().template size<>(),
                2.f * unit<scalar>::mm, tol);

    cl_state.template set_constraint<constraint::e_user>(0.5f *
                                                         unit<scalar>::mm);
    ASSERT_NEAR(cl_state.constraints().template size<>(),
                0.5f * unit<scalar>::mm, tol);

    // Run a few steps
    ASSERT_TRUE(l_stepper.step(propagation));
    // Step constraint to half step size
    ASSERT_TRUE(cl_stepper.step(c_propagation));
    ASSERT_TRUE(cl_stepper.step(c_propagation));

    track = propagation._stepping();
    ASSERT_NEAR(track.pos()[0], constant<scalar>::inv_sqrt2, tol);
    ASSERT_NEAR(track.pos()[1], constant<scalar>::inv_sqrt2, tol);
    ASSERT_NEAR(track.pos()[2], 0.f, tol);

    c_track = c_propagation._stepping();
    ASSERT_NEAR(c_track.pos()[0], constant<scalar>::inv_sqrt2, tol);
    ASSERT_NEAR(c_track.pos()[1], constant<scalar>::inv_sqrt2, tol);
    ASSERT_NEAR(c_track.pos()[2], 0.f, tol);

    ASSERT_TRUE(l_stepper.step(propagation));

    track = propagation._stepping();
    ASSERT_NEAR(track.pos()[0], constant<scalar>::sqrt2, tol);
    ASSERT_NEAR(track.pos()[1], constant<scalar>::sqrt2, tol);
    ASSERT_NEAR(track.pos()[2], 0.f, tol);
}

// This tests the base functionality of the Runge-Kutta stepper
GTEST_TEST(detray_propagator, rk_stepper) {
    using namespace step;

    // RK stepper configurations
    constexpr unsigned int theta_steps = 100u;
    constexpr unsigned int phi_steps = 100u;
    constexpr unsigned int rk_steps = 100u;

    // Constant magnetic field
    vector3 B{1.f * unit<scalar>::T, 1.f * unit<scalar>::T,
              1.f * unit<scalar>::T};
    mag_field_t mag_field(
        typename mag_field_t::backend_t::configuration_t{B[0], B[1], B[2]});

    // RK stepper
    rk_stepper_t rk_stepper;
    crk_stepper_t crk_stepper;

    // Set origin position of tracks
    const point3 ori{0.f, 0.f, 0.f};
    const scalar p_mag{10.f * unit<scalar>::GeV};

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

        crk_state.template set_constraint<constraint::e_user>(0.5f *
                                                              unit<scalar>::mm);
        n_state._step_size = 1.f * unit<scalar>::mm;
        cn_state._step_size = 1.f * unit<scalar>::mm;
        ASSERT_NEAR(crk_state.constraints().template size<>(),
                    0.5f * unit<scalar>::mm, tol);

        for (unsigned int i_s = 0u; i_s < rk_steps; i_s++) {
            rk_stepper.step(propagation);
            crk_stepper.step(c_propagation);
            crk_stepper.step(c_propagation);
        }

        // get relative error by dividing error with path length
        ASSERT_NEAR(rk_state.path_length(), crk_state.path_length(), tol);
        ASSERT_NEAR(getter::norm(rk_state().pos() - crk_state().pos()) /
                        rk_state.path_length(),
                    0.f, tol);

        const auto helix_pos = helix(rk_state.path_length());
        const auto forward_pos = rk_state().pos();
        const point3 forward_relative_error{(1.f / rk_state.path_length()) *
                                            (forward_pos - helix_pos)};

        // Make sure that relative error is smaller than epsion
        EXPECT_NEAR(getter::norm(forward_relative_error), 0.f, tol);

        // Roll the same track back to the origin
        // Use the same path length, since there is no overstepping
        const scalar path_length = rk_state.path_length();
        n_state._step_size *= -unit<scalar>::mm;
        cn_state._step_size *= -unit<scalar>::mm;
        for (unsigned int i_s = 0u; i_s < rk_steps; i_s++) {
            rk_stepper.step(propagation);
            crk_stepper.step(c_propagation);
            crk_stepper.step(c_propagation);
        }

        ASSERT_NEAR(rk_state.path_length(), crk_state.path_length(), tol);
        ASSERT_NEAR(getter::norm(rk_state().pos() - crk_state().pos()) /
                        (2.f * path_length),
                    0.f, tol);

        const auto backward_pos = rk_state().pos();
        const point3 backward_relative_error{1.f / (2.f * path_length) *
                                             (backward_pos - ori)};

        // Make sure that relative error is smaller than epsion
        EXPECT_NEAR(getter::norm(backward_relative_error), 0.f, tol);
    }
}
