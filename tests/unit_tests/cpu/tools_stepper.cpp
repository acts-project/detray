/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// detray include(s)
#include "detray/definitions/bfield_backends.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/tracks.hpp"

// google-test include(s)
#include <gtest/gtest.h>

using namespace detray;
using vector2 = test::vector2;
using vector3 = test::vector3;
using point3 = test::point3;
using transform3 = test::transform3;
using matrix_operator = standard_matrix_operator<scalar>;

/// Runge-Kutta stepper
template <typename bfield_t>
using rk_stepper_t = rk_stepper<typename bfield_t::view_t, transform3>;
template <typename bfield_t>
using crk_stepper_t =
    rk_stepper<typename bfield_t::view_t, transform3, constrained_step<>>;

namespace {

constexpr scalar tol{1e-3f};

template <typename bfield_t>
inline bfield_t load_field() {
    detray::io::detail::file_handle file(
        !std::getenv("DETRAY_BFIELD_FILE") ? ""
                                           : std::getenv("DETRAY_BFIELD_FILE"),
        std::ios_base::binary | std::ios_base::in);

    bfield_t field(*file);

    return field;
}

// dummy navigation struct
struct nav_state {
    scalar operator()() const { return _step_size; }
    inline auto current_object() const -> dindex { return dindex_invalid; }
    inline auto tolerance() const -> scalar { return tol; }

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

    // Constant magnetic field
    using bfield_t = covfie::field<bfield::const_bknd_t>;

    vector3 B{1.f * unit<scalar>::T, 1.f * unit<scalar>::T,
              1.f * unit<scalar>::T};
    bfield_t hom_bfield(covfie::make_parameter_pack(
        typename bfield_t::backend_t::configuration_t{B[0], B[1], B[2]}));

    // RK stepper
    rk_stepper_t<bfield_t> rk_stepper;
    crk_stepper_t<bfield_t> crk_stepper;

    // RK stepper configurations
    constexpr unsigned int rk_steps = 100u;
    constexpr scalar stepsize_constr{0.5f * unit<scalar>::mm};

    // Track generator configuration
    const scalar p_mag{10.f * unit<scalar>::GeV};
    constexpr unsigned int theta_steps = 100u;
    constexpr unsigned int phi_steps = 100u;

    // Iterate through uniformly distributed momentum directions
    for (auto track :
         uniform_track_generator<free_track_parameters<transform3>>(
             phi_steps, theta_steps, p_mag)) {
        // Generate track state used for propagation with constrained step size
        free_track_parameters c_track(track);

        // helix trajectory
        detail::helix helix(track, &B);

        // RK Stepping into forward direction
        prop_state<rk_stepper_t<bfield_t>::state, nav_state> propagation{
            rk_stepper_t<bfield_t>::state{track, hom_bfield}, nav_state{}};
        prop_state<crk_stepper_t<bfield_t>::state, nav_state> c_propagation{
            crk_stepper_t<bfield_t>::state{c_track, hom_bfield}, nav_state{}};

        rk_stepper_t<bfield_t>::state &rk_state = propagation._stepping;
        crk_stepper_t<bfield_t>::state &crk_state = c_propagation._stepping;

        // Retrieve the navigation states
        nav_state &n_state = propagation._navigation;
        nav_state &cn_state = c_propagation._navigation;

        // Set step size constraint to half the nominal step size =>
        // crk_stepper will need twice as many steps
        crk_state.template set_constraint<constraint::e_user>(stepsize_constr);
        ASSERT_NEAR(crk_state.constraints().template size<>(),
                    0.5f * unit<scalar>::mm, tol);

        // Reset step size in the navigation state to a positive value
        n_state._step_size = 1.f * unit<scalar>::mm;
        cn_state._step_size = 1.f * unit<scalar>::mm;
        for (unsigned int i_s = 0u; i_s < rk_steps; i_s++) {
            rk_stepper.step(propagation);
            crk_stepper.step(c_propagation);
            crk_stepper.step(c_propagation);
        }

        // Check that both steppers arrive at the same point
        // Get relative error by dividing error with path length
        ASSERT_NEAR(rk_state.path_length(), crk_state.path_length(), tol);
        ASSERT_NEAR(getter::norm(rk_state().pos() - crk_state().pos()) /
                        rk_state.path_length(),
                    0.f, tol);

        // Check that the stepper position lies on the truth helix
        const auto helix_pos = helix(rk_state.path_length());
        const auto forward_pos = rk_state().pos();
        const point3 forward_relative_error{(1.f / rk_state.path_length()) *
                                            (forward_pos - helix_pos)};

        // Make sure that relative error is smaller than the tolerance
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

        // Should arrive back at track origin, where path length is zero
        ASSERT_NEAR(rk_state.path_length(), 0.f, tol);
        ASSERT_NEAR(crk_state.path_length(), 0.f, tol);

        const point3 backward_relative_error{1.f / (2.f * path_length) *
                                             (rk_state().pos())};
        // Make sure that relative error is smaller than the tolerance
        EXPECT_NEAR(getter::norm(backward_relative_error), 0.f, tol);

        // The constrained stepper should be at the same position now
        ASSERT_NEAR(getter::norm(rk_state().pos() - crk_state().pos()) /
                        (2.f * path_length),
                    0.f, tol);
    }
}

/// This tests the base functionality of the Runge-Kutta stepper in an
/// in-homogeneous magnetic field, read from file
TEST(detray_propagator, rk_stepper_inhomogeneous_bfield) {
    using namespace step;

    using bfield_t = covfie::field<bfield::inhom_bknd_t>;

    // Read the magnetic field map
    bfield_t inhom_bfield(load_field<bfield_t>());

    // RK stepper
    rk_stepper_t<bfield_t> rk_stepper;
    crk_stepper_t<bfield_t> crk_stepper;

    // RK stepper configurations
    constexpr unsigned int rk_steps = 100u;
    constexpr scalar stepsize_constr{0.5f * unit<scalar>::mm};

    // Track generator configuration
    const scalar p_mag{10.f * unit<scalar>::GeV};
    constexpr unsigned int theta_steps = 100u;
    constexpr unsigned int phi_steps = 100u;

    // Iterate through uniformly distributed momentum directions
    for (auto track :
         uniform_track_generator<free_track_parameters<transform3>>(
             phi_steps, theta_steps, p_mag)) {
        // Generate track state used for propagation with constrained step size
        free_track_parameters c_track(track);

        // RK Stepping into forward direction
        prop_state<rk_stepper_t<bfield_t>::state, nav_state> propagation{
            rk_stepper_t<bfield_t>::state{track, inhom_bfield}, nav_state{}};
        prop_state<crk_stepper_t<bfield_t>::state, nav_state> c_propagation{
            crk_stepper_t<bfield_t>::state{c_track, inhom_bfield}, nav_state{}};

        rk_stepper_t<bfield_t>::state &rk_state = propagation._stepping;
        crk_stepper_t<bfield_t>::state &crk_state = c_propagation._stepping;

        // Retrieve the navigation states
        nav_state &n_state = propagation._navigation;
        nav_state &cn_state = c_propagation._navigation;

        crk_state.template set_constraint<constraint::e_user>(stepsize_constr);
        ASSERT_NEAR(crk_state.constraints().template size<>(),
                    0.5f * unit<scalar>::mm, tol);

        n_state._step_size = 1.f * unit<scalar>::mm;
        cn_state._step_size = 1.f * unit<scalar>::mm;
        for (unsigned int i_s = 0u; i_s < rk_steps; i_s++) {
            rk_stepper.step(propagation);
            crk_stepper.step(c_propagation);
            crk_stepper.step(c_propagation);
        }

        // Make sure the steppers moved
        const scalar path_length{rk_state.path_length()};
        ASSERT_TRUE(path_length > 0.f);
        ASSERT_TRUE(crk_state.path_length() > 0.f);

        // Roll the same track back to the origin
        // Use the same path length, since there is no overstepping
        n_state._step_size *= -unit<scalar>::mm;
        cn_state._step_size *= -unit<scalar>::mm;
        for (unsigned int i_s = 0u; i_s < rk_steps; i_s++) {
            rk_stepper.step(propagation);
            crk_stepper.step(c_propagation);
            crk_stepper.step(c_propagation);
        }

        // Should arrive back at track origin, where path length is zero
        ASSERT_NEAR(rk_state.path_length(), 0.f, tol);
        ASSERT_NEAR(crk_state.path_length(), 0.f, tol);

        const point3 backward_relative_error{1.f / (2.f * path_length) *
                                             (rk_state().pos())};
        // Make sure that relative error is smaller than the tolerance
        EXPECT_NEAR(getter::norm(backward_relative_error), 0.f, tol);

        // The constrained stepper should be at the same position now
        ASSERT_NEAR(getter::norm(rk_state().pos() - crk_state().pos()) /
                        (2.f * path_length),
                    0.f, tol);
    }
}
