/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/detectors/detector_metadata.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/inspectors.hpp"

// Covfie include(s).
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>

// Vecmem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// google-test include(s).
#include <gtest/gtest.h>

#include <iostream>

using namespace detray;

namespace env {

// Type Definitions
template <typename mask_shape_t>
using tel_detector_t =
    detector<detector_registry::template telescope_detector<mask_shape_t>,
             covfie::field>;

constexpr scalar tol{1e-2f};

// Surface material
const material<scalar> mat = silicon<scalar>();
const scalar thickness{2.f * unit<scalar>::mm};

// memory resource
vecmem::host_memory_resource resource;

}  // namespace env

/// Test the path correction on a rectangular surface (cartesian coordinates)
TEST(path_correction, cartesian2D) {
    using detector_t = env::tel_detector_t<rectangle2D<>>;
    using mag_field_t = detector_t::bfield_type;

    using matrix_operator = standard_matrix_operator<scalar>;
    using transform3 = typename detector_t::transform3;
    using vector2 = __plugin::vector2<scalar>;
    using vector3 = typename transform3::vector3;

    using navigator_t = navigator<detector_t, navigation::print_inspector>;
    // Re-init the volume at every step, so that surface will be put into the
    // navigation cache
    using crk_stepper_t = rk_stepper<mag_field_t::view_t, transform3,
                                     constrained_step<>, always_init>;
    using actor_chain_t = actor_chain<
        dtuple, target_aborter, pathlimit_aborter, propagation::print_inspector,
        parameter_transporter<transform3>, parameter_resetter<transform3>>;
    using propagator_t = propagator<crk_stepper_t, navigator_t, actor_chain_t>;

    // B field
    const vector3 B{0.f, 0.f, 1.f * unit<scalar>::T};
    mag_field_t mag_field(
        mag_field_t::backend_t::configuration_t{B[0], B[1], B[2]});

    // Generate track starting point
    vector2 local{2.f, 3.f};
    vector3 mom{1.f * unit<scalar>::MeV, 0.f, 0.f};
    scalar time{0.f};
    scalar q{-1.f};

    // Path radius and path length per turn
    const scalar r{getter::perp(mom) / getter::norm(B)};
    const scalar S{2.f * constant<scalar>::pi * r};
    // Maximum distance to propagate
    const scalar path_limit{1.01f * S};

    // Use rectangle surface
    mask<rectangle2D<>> rectangle{0u, 0.9f * r * unit<scalar>::mm,
                                  0.9f * r * unit<scalar>::mm};

    // Used to build the telescope detector along x-axis
    detail::ray<transform3> pilot_traj{{0.f, 0.f, 0.f}, time, mom, q};

    // Build telescope detector with a single rectangle surfaces and an envelope
    // between the surface and the portals that allows to propagate a full turn
    const auto det = create_telescope_detector(
        env::resource, std::move(mag_field), rectangle, 1u, S, env::mat,
        env::thickness, pilot_traj, 1.75f * r);

    // bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0u) = local[0];
    getter::element(bound_vector, e_bound_loc1, 0u) = local[1];
    getter::element(bound_vector, e_bound_phi, 0u) = getter::phi(mom);
    getter::element(bound_vector, e_bound_theta, 0u) = getter::theta(mom);
    getter::element(bound_vector, e_bound_qoverp, 0u) = q / getter::norm(mom);
    getter::element(bound_vector, e_bound_time, 0u) = time;

    // bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov, e_bound_loc0, e_bound_loc0) = 1.f;
    getter::element(bound_cov, e_bound_loc1, e_bound_loc1) = 1.f;
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 1.f;
    // Note: Set theta error as ZERO, to constrain the loc1 divergence
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.f;
    getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 1.f;
    getter::element(bound_cov, e_bound_time, e_bound_time) = 1.f;

    // bound track parameter
    geometry::barcode bcd{};
    bcd.set_volume(0u).set_id(surface_id::e_sensitive).set_index(0u);
    const bound_track_parameters<transform3> bound_param0(bcd, bound_vector,
                                                          bound_cov);

    // Actors
    target_aborter::state targeter{bcd};
    pathlimit_aborter::state pathlimit_state{path_limit};
    propagation::print_inspector::state print_insp_state{};
    parameter_transporter<transform3>::state bound_updater{};
    parameter_resetter<transform3>::state rst{};

    propagator_t p({}, {});
    propagator_t::state propagation(bound_param0, det.get_bfield(), det);
    crk_stepper_t::state &crk_state = propagation._stepping;

    // The track is strongly bent, allow for overstepping
    crk_state().set_overstep_tolerance(-1.f * unit<scalar>::um);

    // Retrieve navigation information
    auto &debug_printer = propagation._navigation.inspector();

    // Run propagator
    ASSERT_TRUE(p.propagate(
        propagation, std::tie(targeter, pathlimit_state, print_insp_state,
                              bound_updater, rst)));

    // Make sure, the track was propagated to the end
    ASSERT_NEAR(crk_state.path_length(), S, env::tol)
        << print_insp_state.to_string() << debug_printer.to_string();
    EXPECT_NEAR(crk_state().pos()[0], 0.f, env::tol);
    EXPECT_NEAR(crk_state().pos()[1], 2.f, env::tol);
    EXPECT_NEAR(crk_state().pos()[2], 3.f, env::tol);

    // Bound state after one turn propagation
    const auto bound_param1 = crk_state._bound_params;
    const auto bound_vec0 = bound_param0.vector();
    const auto bound_vec1 = bound_param1.vector();
    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check if the bound state stays the same after one turn propagation

    // vector
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec0, i, 0u),
                    matrix_operator().element(bound_vec1, i, 0u), env::tol);
    }

    // covariance
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j), env::tol);
        }
    }

    // Print the matrix elements
    /*for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            printf("%f ", matrix_operator().element(bound_cov0, i, j));
        }
        printf("\n");
    }
    printf("\n");

    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            printf("%f ", matrix_operator().element(bound_cov1, i, j));
        }
        printf("\n");
    }
    printf("\n");*/
}

/// Test the path correction on a ring surface (polar coordinates)
TEST(path_correction, polar) {

    using detector_t = env::tel_detector_t<ring2D<>>;
    using mag_field_t = detector_t::bfield_type;

    using matrix_operator = standard_matrix_operator<scalar>;
    using transform3 = typename detector_t::transform3;
    using vector2 = __plugin::vector2<scalar>;
    using vector3 = typename transform3::vector3;

    using navigator_t = navigator<detector_t, navigation::print_inspector>;
    // Re-init the volume at every step, so that surface will be put into the
    // navigation cache
    using crk_stepper_t = rk_stepper<mag_field_t::view_t, transform3,
                                     constrained_step<>, always_init>;
    using actor_chain_t = actor_chain<
        dtuple, target_aborter, pathlimit_aborter, propagation::print_inspector,
        parameter_transporter<transform3>, parameter_resetter<transform3>>;
    using propagator_t = propagator<crk_stepper_t, navigator_t, actor_chain_t>;

    // B field
    const vector3 B{0.f, 0.f, 1.f * unit<scalar>::T};
    mag_field_t mag_field(
        mag_field_t::backend_t::configuration_t{B[0], B[1], B[2]});

    // Generate track starting point
    vector2 local{2.f, constant<scalar>::pi / 6.f};
    vector3 mom{1.f * unit<scalar>::MeV, 0.f, 0.f};
    scalar time{0.f};
    scalar q{-1.f};

    // Path radius and path length per turn
    const scalar r{getter::perp(mom) / getter::norm(B)};
    const scalar S{2.f * constant<scalar>::pi * r};
    // Maximum distance to propagate
    const scalar path_limit{1.01f * S};

    // Use ring surface
    const scalar r_low{0.f * unit<scalar>::mm};
    const scalar r_high{0.9f * r * unit<scalar>::mm};
    mask<ring2D<>> ring{0u, r_low, r_high};

    // Used to build the telescope detector along x-axis
    detail::ray<transform3> pilot_traj{{0.f, 0.f, 0.f}, time, mom, q};

    // Build telescope detector with a single ring surfaces and an envelope
    // between the surface and the portals that allows to propagate a full turn
    const auto det = create_telescope_detector(
        env::resource, std::move(mag_field), ring, 1u, S, env::mat,
        env::thickness, pilot_traj, 1.75f * r);

    // bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0u) = local[0];
    getter::element(bound_vector, e_bound_loc1, 0u) = local[1];
    getter::element(bound_vector, e_bound_phi, 0u) = getter::phi(mom);
    getter::element(bound_vector, e_bound_theta, 0u) = getter::theta(mom);
    getter::element(bound_vector, e_bound_qoverp, 0u) = q / getter::norm(mom);
    getter::element(bound_vector, e_bound_time, 0u) = time;

    // bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov, e_bound_loc0, e_bound_loc0) = 1.f;
    getter::element(bound_cov, e_bound_loc1, e_bound_loc1) = 1.f;
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 1.f;
    // Note: Set theta error as ZERO, to constrain the loc1 divergence
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.f;
    getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 1.f;
    getter::element(bound_cov, e_bound_time, e_bound_time) = 1.f;

    // bound track parameter
    geometry::barcode bcd{};
    bcd.set_volume(0u).set_id(surface_id::e_sensitive).set_index(0u);
    const bound_track_parameters<transform3> bound_param0(bcd, bound_vector,
                                                          bound_cov);

    // Actors
    target_aborter::state targeter{bcd};
    pathlimit_aborter::state pathlimit_state{path_limit};
    propagation::print_inspector::state print_insp_state{};
    parameter_transporter<transform3>::state bound_updater{};
    parameter_resetter<transform3>::state rst{};

    propagator_t p({}, {});
    propagator_t::state propagation(bound_param0, det.get_bfield(), det);
    crk_stepper_t::state &crk_state = propagation._stepping;

    // The track is strongly bent, allow for overstepping
    crk_state().set_overstep_tolerance(-1.f * unit<scalar>::um);

    // Retrieve navigation information
    auto &debug_printer = propagation._navigation.inspector();

    // Run propagator
    ASSERT_TRUE(p.propagate(
        propagation, std::tie(targeter, pathlimit_state, print_insp_state,
                              bound_updater, rst)));

    // Make sure, the track was propagated to the end
    ASSERT_NEAR(crk_state.path_length(), S, env::tol)
        << print_insp_state.to_string() << debug_printer.to_string();
    EXPECT_NEAR(crk_state().pos()[0], 0.f, env::tol);
    EXPECT_NEAR(crk_state().pos()[1],
                2.f * std::cos(constant<scalar>::pi / 6.f), env::tol);
    EXPECT_NEAR(crk_state().pos()[2],
                2.f * std::sin(constant<scalar>::pi / 6.f), env::tol);

    // Bound state after one turn propagation
    const auto bound_param1 = crk_state._bound_params;
    const auto bound_vec0 = bound_param0.vector();
    const auto bound_vec1 = bound_param1.vector();
    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check if the bound state stays the same after one turn propagation

    // vector
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec0, i, 0u),
                    matrix_operator().element(bound_vec1, i, 0u), env::tol);
    }

    // covariance
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j), env::tol);
        }
    }
}

/// Test the path correction on a cylindrical surface (cylinder coordinates)
TEST(path_correction, cylindrical) {

    using detector_t = env::tel_detector_t<cylinder2D<>>;
    using mag_field_t = detector_t::bfield_type;

    using matrix_operator = standard_matrix_operator<scalar>;
    using transform3 = typename detector_t::transform3;
    using vector2 = __plugin::vector2<scalar>;
    using vector3 = typename transform3::vector3;

    using navigator_t = navigator<detector_t, navigation::print_inspector>;
    // Re-init the volume at every step, so that surface will be put into the
    // navigation cache
    using crk_stepper_t = rk_stepper<mag_field_t::view_t, transform3,
                                     constrained_step<>, always_init>;
    using actor_chain_t = actor_chain<
        dtuple, target_aborter, pathlimit_aborter, propagation::print_inspector,
        parameter_transporter<transform3>, parameter_resetter<transform3>>;
    using propagator_t = propagator<crk_stepper_t, navigator_t, actor_chain_t>;

    // B field
    const vector3 B{0.f, 0.f, 1.f * unit<scalar>::T};
    mag_field_t mag_field(
        mag_field_t::backend_t::configuration_t{B[0], B[1], B[2]});

    // Generate track starting point: loc phi=0 with direction pointing upwards
    // along global y-axis and x- and y-axis get switched
    vector2 local{0.f, 0.f};
    vector3 mom{0.f, 1.f * unit<scalar>::MeV, 0.f};
    scalar time{0.f};
    scalar q{-1.f};

    // Path radius and path length per turn
    const scalar r{getter::perp(mom) / getter::norm(B)};
    const scalar S{2.f * constant<scalar>::pi * r};
    // Maximum distance to propagate
    const scalar path_limit{1.01f * S};

    // Use 2-dimensional cylinder surface
    const scalar max_r{0.5f * r};
    const scalar n_half_z{-0.9f * r};
    const scalar p_half_z{0.9f * r};
    mask<cylinder2D<>> cyl{0u, max_r, n_half_z, p_half_z};

    // Used to build the telescope detector along x-axis, this rotates the
    // cylinder so that its axis lies along x
    detail::ray<transform3> pilot_traj{
        {0.f, 0.f, 0.f}, time, {1.f, 0.f, 0.f}, q};

    // Build telescope detector with a single cylinder surfaces and an envelope
    // between the surface and the portals that allows to propagate a full turn
    const auto det = create_telescope_detector(
        env::resource, std::move(mag_field), cyl, 1u, S, env::mat,
        env::thickness, pilot_traj, 1.75f * r + max_r);

    // bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0u) = local[0];
    getter::element(bound_vector, e_bound_loc1, 0u) = local[1];
    getter::element(bound_vector, e_bound_phi, 0u) = getter::phi(mom);
    getter::element(bound_vector, e_bound_theta, 0u) = getter::theta(mom);
    getter::element(bound_vector, e_bound_qoverp, 0u) = q / getter::norm(mom);
    getter::element(bound_vector, e_bound_time, 0u) = time;

    // bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov, e_bound_loc0, e_bound_loc0) = 1.f;
    getter::element(bound_cov, e_bound_loc1, e_bound_loc1) = 1.f;
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 1.f;
    // Note: Set theta error as ZERO, to constrain the loc1 divergence
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.f;
    getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 1.f;
    getter::element(bound_cov, e_bound_time, e_bound_time) = 1.f;

    // bound track parameter
    geometry::barcode bcd{};
    bcd.set_volume(0u).set_id(surface_id::e_sensitive).set_index(0u);
    const bound_track_parameters<transform3> bound_param0(bcd, bound_vector,
                                                          bound_cov);

    // Actors
    target_aborter::state targeter{bcd};
    pathlimit_aborter::state pathlimit_state{path_limit};
    propagation::print_inspector::state print_insp_state{};
    parameter_transporter<transform3>::state bound_updater{};
    parameter_resetter<transform3>::state rst{};

    propagator_t p({}, {});
    propagator_t::state propagation(bound_param0, det.get_bfield(), det);
    crk_stepper_t::state &crk_state = propagation._stepping;

    // The track is strongly bent, allow for overstepping
    crk_state().set_overstep_tolerance(-1.f * unit<scalar>::um);

    // RK stepper and its state
    EXPECT_NEAR(crk_state().pos()[0], 0.f, env::tol);    // x
    EXPECT_NEAR(crk_state().pos()[1], max_r, env::tol);  // y
    EXPECT_NEAR(crk_state().pos()[2], 0.f, env::tol);    // z

    // Retrieve navigation information
    auto &debug_printer = propagation._navigation.inspector();

    // Run propagator
    ASSERT_TRUE(p.propagate(
        propagation, std::tie(targeter, pathlimit_state, print_insp_state,
                              bound_updater, rst)))
        << print_insp_state.to_string() << debug_printer.to_string();

    // Make sure, the track was propagated to the end
    ASSERT_NEAR(crk_state.path_length(), S, env::tol)
        << print_insp_state.to_string() << debug_printer.to_string();
    EXPECT_NEAR(crk_state().pos()[0], 0.f, env::tol);    // x
    EXPECT_NEAR(crk_state().pos()[1], max_r, env::tol);  // y
    EXPECT_NEAR(crk_state().pos()[2], 0.f, env::tol);    // z

    // Bound state after one turn propagation
    const auto bound_param1 = crk_state._bound_params;
    const auto bound_vec0 = bound_param0.vector();
    const auto bound_vec1 = bound_param1.vector();
    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check if the bound state stays the same after one turn propagation

    // vector
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec0, i, 0u),
                    matrix_operator().element(bound_vec1, i, 0u), env::tol);
    }

    // covariance
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j), env::tol);
        }
    }
}