/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/geometry/barcode.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s)
#include "detray/test/utils/detectors/build_telescope_detector.hpp"
#include "detray/test/utils/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// google-test include(s).
#include <gtest/gtest.h>

using namespace detray;

// Algebra types
using test_algebra = test::algebra;
using scalar = test::scalar;
using point2 = test::point2;
using vector3 = test::vector3;

// Test class for the backward propagation
// Input tuple: < std::vector<plane positions>, tolerance >
class BackwardPropagation
    : public ::testing::TestWithParam<std::tuple<std::vector<scalar>, scalar>> {
};

TEST_P(BackwardPropagation, backward_propagation) {

    const scalar tol = std::get<1>(GetParam());

    vecmem::host_memory_resource host_mr;

    // Build in x-direction from given module positions
    detail::ray<test_algebra> traj{{0.f, 0.f, 0.f}, 0.f, {1.f, 0.f, 0.f}, -1.f};
    std::vector<scalar> positions = std::get<0>(GetParam());

    tel_det_config<test_algebra, rectangle2D> tel_cfg{200.f * unit<scalar>::mm,
                                                      200.f * unit<scalar>::mm};
    tel_cfg.positions(positions).pilot_track(traj).mat_thickness(
        10.f * unit<scalar>::mm);

    // Build telescope detector with rectangular planes
    const auto [det, names] = build_telescope_detector(host_mr, tel_cfg);

    // Create b field
    using bfield_t = bfield::const_field_t<scalar>;
    vector3 B{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
              1.f * unit<scalar>::T};
    const bfield_t hom_bfield = bfield::create_const_field<scalar>(B);

    using navigator_t = navigator<decltype(det)>;
    using rk_stepper_t = rk_stepper<bfield_t::view_t, test_algebra>;
    using actor_chain_t =
        actor_chain<parameter_transporter<test_algebra>,
                    pointwise_material_interactor<test_algebra>,
                    parameter_resetter<test_algebra>>;
    using propagator_t = propagator<rk_stepper_t, navigator_t, actor_chain_t>;

    // Particle hypothesis
    pdg_particle<scalar> ptc = muon<scalar>();

    // Bound vector
    bound_parameters_vector<test_algebra> bound_vector{};
    bound_vector.set_theta(constant<scalar>::pi_2);
    bound_vector.set_qop(ptc.charge() / (1.f * unit<scalar>::GeV));

    // Bound covariance
    auto bound_cov = matrix::identity<
        typename bound_track_parameters<test_algebra>::covariance_type>();

    // Bound track parameter
    const bound_track_parameters<test_algebra> bound_param0(
        geometry::barcode{}.set_index(0u), bound_vector, bound_cov);

    // Actors
    parameter_transporter<test_algebra>::state bound_updater{};
    pointwise_material_interactor<test_algebra>::state interactor{};
    parameter_resetter<test_algebra>::state rst{};

    propagation::config prop_cfg{};
    prop_cfg.stepping.rk_error_tol = 1e-12f * unit<float>::mm;
    prop_cfg.navigation.overstep_tolerance = -100.f * unit<float>::um;
    propagator_t p{prop_cfg};

    // Forward state
    propagator_t::state fw_state(bound_param0, hom_bfield, det,
                                 prop_cfg.context);
    fw_state.set_particle(ptc);
    fw_state.do_debug = true;

    // Run propagator
    p.propagate(fw_state, detray::tie(bound_updater, interactor, rst));

    // Print the debug stream
    // std::cout << fw_state.debug_stream.str() << std::endl;

    // Bound state after propagation
    const auto& bound_param1 = fw_state._stepping.bound_params();

    // Check if the track reaches the final surface
    EXPECT_EQ(bound_param0.surface_link().volume(), 4095u);
    EXPECT_EQ(bound_param0.surface_link().index(), 0u);
    EXPECT_EQ(bound_param1.surface_link().volume(), 0u);
    EXPECT_EQ(bound_param1.surface_link().id(), surface_id::e_sensitive);
    EXPECT_EQ(bound_param1.surface_link().index(), positions.size() - 1u);

    // Backward state
    propagator_t::state bw_state(bound_param1, hom_bfield, det,
                                 prop_cfg.context);
    bw_state.set_particle(ptc);
    bw_state.do_debug = true;
    bw_state._navigation.set_direction(navigation::direction::e_backward);

    // Run propagator
    p.propagate(bw_state, detray::tie(bound_updater, interactor, rst));

    // Print the debug stream
    // std::cout << bw_state.debug_stream.str() << std::endl;

    // Bound state after propagation
    const auto& bound_param2 = bw_state._stepping.bound_params();

    // Check if the track reaches the initial surface
    EXPECT_EQ(bound_param2.surface_link().volume(), 0u);
    EXPECT_EQ(bound_param2.surface_link().id(), surface_id::e_sensitive);
    EXPECT_EQ(bound_param2.surface_link().index(), 0u);

    const auto bound_vec0 = bound_param0.vector();
    const auto bound_vec2 = bound_param2.vector();

    // Check vector
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        EXPECT_NEAR(getter::element(bound_vec0, i, 0),
                    getter::element(bound_vec2, i, 0), tol);
    }

    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Some sanity checks
    EXPECT_TRUE(bound_param0.p(ptc.charge()) > bound_param1.p(ptc.charge()));
    EXPECT_TRUE(bound_param2.p(ptc.charge()) > bound_param1.p(ptc.charge()));

    EXPECT_TRUE(getter::element(bound_cov1, e_bound_qoverp, e_bound_qoverp) >
                getter::element(bound_cov0, e_bound_qoverp, e_bound_qoverp));
    EXPECT_TRUE(getter::element(bound_cov1, e_bound_theta, e_bound_theta) >
                getter::element(bound_cov0, e_bound_theta, e_bound_theta));
    EXPECT_TRUE(getter::element(bound_cov1, e_bound_phi, e_bound_phi) >
                getter::element(bound_cov0, e_bound_phi, e_bound_phi));
}

INSTANTIATE_TEST_SUITE_P(
    telescope, BackwardPropagation,
    ::testing::Values(
        std::make_tuple(std::vector<scalar>{0.f}, 1e-5f),
        std::make_tuple(std::vector<scalar>{0.f, 10.f}, 1e-3f),
        std::make_tuple(std::vector<scalar>{0.f, 10.f, 20.f, 30.f, 40.f, 50.f,
                                            60.f, 70.f, 80.f, 90.f, 100.f},
                        1e-2f)));
