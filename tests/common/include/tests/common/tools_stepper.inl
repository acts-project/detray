/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// detray include(s)
#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/propagator/track.hpp"
#include "tests/common/tools/helix_gun.hpp"

// google-test include(s)
#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;
using namespace __plugin;

constexpr scalar epsilon = 1e-4;
constexpr scalar path_limit = 100 * unit_constants::cm;

// dummy navigation struct
struct nav_state {
    scalar operator()() const { return _step_size; }
    inline void set_full_trust() {}
    inline void set_high_trust() {}
    inline void set_fair_trust() {}
    inline void set_no_trust() {}
    inline bool abort() { return false; }

    scalar _step_size = 1. * unit_constants::mm;
};
nav_state n_state{};

// This tests the base functionality of the line stepper
TEST(ALGEBRA_PLUGIN, line_stepper) {
    using namespace step;

    // Line stepper with and without constrained stepping
    using line_stepper_t = line_stepper<free_track_parameters>;
    using cline_stepper_t =
        line_stepper<free_track_parameters, constrained_step<>>;

    point3<scalar> pos{0., 0., 0.};
    vector3<scalar> mom{1., 1., 0.};
    free_track_parameters traj(pos, 0, mom, -1);
    free_track_parameters c_traj(pos, 0, mom, -1);

    line_stepper_t l_stepper;
    cline_stepper_t cl_stepper;
    line_stepper_t::state l_state(traj);
    cline_stepper_t::state cl_state(c_traj);

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
    ASSERT_TRUE(l_stepper.step(l_state, n_state));
    // Step constraint to half step size
    ASSERT_TRUE(cl_stepper.step(cl_state, n_state));
    ASSERT_TRUE(cl_stepper.step(cl_state, n_state));

    ASSERT_FLOAT_EQ(traj.pos()[0], 1. / std::sqrt(2));
    ASSERT_FLOAT_EQ(traj.pos()[1], 1. / std::sqrt(2));
    ASSERT_FLOAT_EQ(traj.pos()[2], 0.);

    ASSERT_FLOAT_EQ(c_traj.pos()[0], 1. / std::sqrt(2));
    ASSERT_FLOAT_EQ(c_traj.pos()[1], 1. / std::sqrt(2));
    ASSERT_FLOAT_EQ(c_traj.pos()[2], 0.);

    // Only step 1.5 mm further, then break
    l_state.set_path_limit(1.5);
    ASSERT_TRUE(l_stepper.step(l_state, n_state));

    ASSERT_FLOAT_EQ(traj.pos()[0], std::sqrt(2));
    ASSERT_FLOAT_EQ(traj.pos()[1], std::sqrt(2));
    ASSERT_FLOAT_EQ(traj.pos()[2], 0.);

    // breaks
    ASSERT_FALSE(l_stepper.step(l_state, n_state));
}

// This tests the base functionality of the Runge-Kutta stepper
TEST(ALGEBRA_PLUGIN, rk_stepper) {
    using namespace step;

    // type definitions
    using vector3 = vector3<scalar>;
    using point3 = point3<scalar>;

    using mag_field_t = constant_magnetic_field<>;
    using rk_stepper_t = rk_stepper<mag_field_t, free_track_parameters>;
    using crk_stepper_t =
        rk_stepper<mag_field_t, free_track_parameters, constrained_step<>>;

    // RK stepper configurations
    constexpr unsigned int theta_steps = 2;
    constexpr unsigned int phi_steps = 2;
    constexpr unsigned int rk_steps = 2;

    // Constant magnetic field
    vector3 B{1 * unit_constants::T, 1 * unit_constants::T,
              1 * unit_constants::T};
    mag_field_t mag_field(B);

    // RK stepper
    rk_stepper_t rk_stepper(mag_field);
    crk_stepper_t crk_stepper(mag_field);

    // Set origin position of tracks
    const point3 ori{0., 0., 0.};
    scalar p_mag = 10;

    // Loops of theta values ]0,pi[
    for (unsigned int itheta = 0; itheta < theta_steps; ++itheta) {
        scalar theta = 0.1 + itheta * (M_PI - 0.1) / theta_steps;
        scalar sin_theta = std::sin(theta);
        scalar cos_theta = std::cos(theta);

        // Loops of phi values [-pi, pi]
        for (unsigned int iphi = 0; iphi < phi_steps; ++iphi) {

            // The direction
            scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
            scalar sin_phi = std::sin(phi);
            scalar cos_phi = std::cos(phi);
            const vector3 dir{cos_phi * sin_theta, sin_phi * sin_theta,
                              cos_theta};

            vector3 mom = p_mag * dir;
            free_track_parameters traj(ori, 0, mom, -1);
            free_track_parameters c_traj(ori, 0, mom, -1);

            // helix gun
            helix_gun helix(traj, &B);

            // RK Stepping into forward direction
            rk_stepper_t::state rk_state(traj);
            crk_stepper_t::state crk_state(c_traj);

            crk_state.template set_constraint<constraint::e_user>(
                0.5 * unit_constants::mm);
            n_state._step_size = 1. * unit_constants::mm;
            ASSERT_NEAR(crk_state.constraints().template size<>(),
                        0.5 * unit_constants::mm, epsilon);

            rk_state.set_path_limit(path_limit);
            crk_state.set_path_limit(path_limit);
            for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
                rk_stepper.step(rk_state, n_state);
                crk_stepper.step(crk_state, n_state);
                crk_stepper.step(crk_state, n_state);
            }

            // get relative error by dividing error with path length
            ASSERT_NEAR(rk_state.dist_to_path_limit(),
                        crk_state.dist_to_path_limit(), epsilon);
            ASSERT_NEAR(getter::norm(rk_state().pos() - crk_state().pos()), 0,
                        epsilon);

            scalar path_accumulated =
                path_limit - rk_state.dist_to_path_limit();
            point3 helix_pos = helix(path_accumulated);
            point3 forward_pos = rk_state().pos();
            vector3 forward_relative_error =
                static_cast<scalar>(1. / path_accumulated) *
                (forward_pos - helix_pos);

            // Make sure that relative error is smaller than epsion
            EXPECT_NEAR(getter::norm(forward_relative_error), 0, epsilon);

            // Roll the same track back to the origin
            n_state._step_size = -1. * unit_constants::mm;
            scalar bw_startig_point = rk_state.dist_to_path_limit();
            for (unsigned int i_s = 0; i_s < rk_steps; i_s++) {
                rk_stepper.step(rk_state, n_state);
                crk_stepper.step(crk_state, n_state);
                crk_stepper.step(crk_state, n_state);
            }

            ASSERT_NEAR(rk_state.dist_to_path_limit(),
                        crk_state.dist_to_path_limit(), epsilon);
            ASSERT_NEAR(getter::norm(rk_state().pos() - crk_state().pos()), 0,
                        epsilon);

            path_accumulated = rk_state.dist_to_path_limit() - bw_startig_point;
            point3 backward_pos = rk_state().pos();
            vector3 backward_relative_error =
                static_cast<scalar>(1. / path_accumulated) *
                (backward_pos - ori);

            // Make sure that relative error is smaller than epsion
            EXPECT_NEAR(getter::norm(backward_relative_error), 0, epsilon);
        }
    }
}

// This tests the covariance transport in rk stepper
TEST(ALGEBRA_PLUGIN, covariance_transport) {

    using vector3 = vector3<scalar>;
    using point3 = point3<scalar>;
    using matrix_operator = standard_matrix_operator<scalar>;
    using mag_field_t = constant_magnetic_field<>;
    using crk_stepper_t =
        rk_stepper<mag_field_t, free_track_parameters, constrained_step<>>;

    // Generate track starting point
    point3 pos{0., 0., 0.};
    scalar time = 0.;
    vector3 mom{0.01, 0.01, 0.01};
    scalar q = -1.;

    free_track_parameters vertex(pos, time, mom, q);

    // Fill covariance matrix
    auto covariance0 =
        matrix_operator().template zero<e_free_size, e_free_size>();

    matrix_operator().element(covariance0, e_free_pos0, e_free_pos0) = 0.1;
    matrix_operator().element(covariance0, e_free_pos1, e_free_pos1) = 0.1;
    matrix_operator().element(covariance0, e_free_pos2, e_free_pos2) = 0.1;
    matrix_operator().element(covariance0, e_free_dir0, e_free_dir0) = 0.1;
    matrix_operator().element(covariance0, e_free_dir1, e_free_dir1) = 0.1;
    matrix_operator().element(covariance0, e_free_dir2, e_free_dir2) = 0.;
    matrix_operator().element(covariance0, e_free_qoverp, e_free_qoverp) = 0.;

    vertex.set_covariance(covariance0);

    // magnetic field
    vector3 B{0, 0, 10. * unit_constants::T};

    // helix gun
    helix_gun helix(vertex, &B);

    // Path length per turn
    scalar S = 2. * getter::norm(mom) / getter::norm(B) * M_PI;

    // RK stepper and its state
    mag_field_t mag_field(B);
    crk_stepper_t::state crk_state(vertex);

    // Decrease tolerance down to 1e-8
    crk_state._tolerance = 1e-8;

    crk_stepper_t crk_stepper(mag_field);

    // Run stepper for one turn
    unsigned int max_steps = 1e4;
    for (unsigned int i = 0; i < max_steps; i++) {

        crk_state.set_constraint(S - crk_state.path_length());

        n_state._step_size = std::numeric_limits<scalar>::max();

        crk_stepper.step(crk_state, n_state);

        if (std::abs(S - crk_state.path_length()) < 1e-6) {
            break;
        }

        // Make sure that we didn't reach the end of for loop
        ASSERT_TRUE(i < max_steps - 1);
    }

    auto jac_transport = crk_state._jac_transport;

    ///// Jacobian check with Helix Gun
    auto J = helix.jacobian(crk_state.path_length());

    for (size_type i = 0; i < e_free_size; i++) {
        for (size_type j = 0; j < e_free_size; j++) {
            EXPECT_NEAR(matrix_operator().element(jac_transport, i, j),
                        matrix_operator().element(J, i, j),
                        crk_state.path_length() * epsilon);
        }
    }

    // Check if covariance doesn't change after one turn
    free_matrix covariance1 = jac_transport * crk_state().covariance() *
                              matrix_operator().transpose(jac_transport);

    for (size_type i = 0; i < e_free_size; i++) {
        for (size_type j = 0; j < e_free_size; j++) {
            EXPECT_NEAR(matrix_operator().element(covariance0, i, j),
                        matrix_operator().element(covariance1, i, j),
                        crk_state.path_length() * epsilon);
        }
    }
}