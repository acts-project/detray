/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/propagator/line_stepper.hpp"

#include "detray/definitions/bfield_backends.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/unbounded.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/test/types.hpp"

// Vecmem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// google-test include(s).
#include <gtest/gtest.h>

using namespace detray;
using matrix_operator = standard_matrix_operator<scalar>;
using transform3 = test::transform3;

constexpr scalar tol{1e-6f};

GTEST_TEST(detray_propagator, covariance_transport) {

    vecmem::host_memory_resource host_mr;

    // Build in x-direction from given module positions
    detail::ray<transform3> traj{{0.f, 0.f, 0.f}, 0.f, {1.f, 0.f, 0.f}, -1.f};
    std::vector<scalar> positions = {0.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f};

    tel_det_config<rectangle2D<>> tel_cfg{200.f * unit<scalar>::mm,
                                          200.f * unit<scalar>::mm};
    tel_cfg.positions(positions).pilot_track(traj);

    // Build telescope detector with unbounded planes
    const auto [det, names] = create_telescope_detector(host_mr, tel_cfg);

    using navigator_t = navigator<decltype(det)>;
    using cline_stepper_t = line_stepper<transform3>;
    using actor_chain_t = actor_chain<dtuple, parameter_transporter<transform3>,
                                      parameter_resetter<transform3>>;
    using propagator_t =
        propagator<cline_stepper_t, navigator_t, actor_chain_t>;

    // Bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0u) = 0.f;
    getter::element(bound_vector, e_bound_loc1, 0u) = 0.f;
    getter::element(bound_vector, e_bound_phi, 0u) = 0.f;
    getter::element(bound_vector, e_bound_theta, 0u) = constant<scalar>::pi_4;
    getter::element(bound_vector, e_bound_qoverp, 0u) = -0.1f;
    getter::element(bound_vector, e_bound_time, 0u) = 0.f;

    // Bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template identity<e_bound_size, e_bound_size>();

    // Note: Set angle error as ZERO, to constrain the loc0 and loc1 divergence
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 0.f;
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.f;

    // Bound track parameter
    const bound_track_parameters<transform3> bound_param0(
        geometry::barcode{}.set_index(0u), bound_vector, bound_cov);

    // Actors
    parameter_transporter<transform3>::state bound_updater{};
    parameter_resetter<transform3>::state rst{};

    propagator_t p({}, {});
    propagator_t::state propagation(bound_param0, det);

    // Run propagator
    p.propagate(propagation, std::tie(bound_updater, rst));

    // Bound state after one turn propagation
    const auto &bound_param1 = propagation._stepping._bound_params;

    // const auto bound_vec0 = bound_param0.vector();
    // const auto bound_vec1 = bound_param1.vector();

    // Check if the track reaches the final surface
    EXPECT_EQ(bound_param0.surface_link().volume(), 4095u);
    EXPECT_EQ(bound_param0.surface_link().index(), 0u);
    EXPECT_EQ(bound_param1.surface_link().volume(), 0u);
    EXPECT_EQ(bound_param1.surface_link().id(), surface_id::e_sensitive);
    EXPECT_EQ(bound_param1.surface_link().index(), 6u);

    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check covaraince
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j), tol);
        }
    }

    // Check covaraince
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j), tol);
        }
    }
}
