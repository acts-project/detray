/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/unbounded.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/axis_rotation.hpp"
#include "tests/common/tools/intersectors/helix_plane_intersector.hpp"

// google-test include(s).
#include <gtest/gtest.h>

using namespace detray;

using matrix_operator = standard_matrix_operator<scalar>;
using transform3 = __plugin::transform3<scalar>;
using vector3 = typename transform3::vector3;

// Setups
constexpr const scalar tolerance = scalar(2e-2);
const vector3 y_axis{0.f, 1.f, 0.f};
const vector3 z_axis{0.f, 0.f, 1.f};
const vector3 B{0.f, 0.f, 1.f * unit<scalar>::T};

// Algebra objects
constexpr const cartesian2<transform3> c2;
constexpr const detail::helix_plane_intersector<transform3> hpi;

/// Test the path correction on a rectangular surface (cartesian coordinates)
TEST(helix_covariance_transport, cartesian2D) {

    // Rectangle surface with enough size
    const detray::mask<detray::unbounded<detray::rectangle2D<>>> rectangle{
        0u, 1e5f * unit<scalar>::mm, 1e5f * unit<scalar>::mm};

    // @NOTE: The test with high energy (>1 GeV) and SINGLE precision fails due
    // to the numerical instability
    free_track_parameters<transform3> free_trk(
        {0.f, 0.f, 0.f}, 0.f, {0.1f * unit<scalar>::GeV, 0.f, 0.f}, -1.f);

    // Path length per turn
    const scalar r{free_trk.p() / getter::norm(B)};
    const scalar S{2.f * constant<scalar>::pi * r};

    // Number of planes and step size between two planes
    const std::size_t n_planes = 10;
    const scalar step_size = S / n_planes;

    // Reference Helix
    detail::helix<transform3> reference_helix(free_trk, &B);

    // Rotation axis in z direction with 20 degree (Can be an arbitrary angle)
    axis_rotation<transform3> axis_rot(z_axis, constant<scalar>::pi / 9.f);

    // Prepare transform matrices
    std::vector<transform3> trfs;
    for (std::size_t i = 0; i < n_planes; i++) {

        const scalar s = step_size * scalar(i);

        // Translation of the new surface
        const auto T = reference_helix(s);

        // Normal vector of the new surface
        auto w = reference_helix.dir(s);
        w = axis_rot(w);

        // Vector on the surface
        const auto v = vector::cross(z_axis, w);

        // Add transform matrix
        trfs.emplace_back(T, w, v);
    }

    // Set the initial bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov0 =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov0, e_bound_loc0, e_bound_loc0) = 1.f;
    getter::element(bound_cov0, e_bound_loc1, e_bound_loc1) = 1.f;
    getter::element(bound_cov0, e_bound_phi, e_bound_phi) = 1.f;
    // Set theta error to zero, to suppress the loc1 (z) divergence
    getter::element(bound_cov0, e_bound_theta, e_bound_theta) = 0.f;
    getter::element(bound_cov0, e_bound_qoverp, e_bound_qoverp) = 1.f;
    getter::element(bound_cov0, e_bound_time, e_bound_time) = 0.f;

    // Copy the covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        bound_cov0;

    // Total path length, just for testing purpose
    scalar total_path_length = 0.f;

    // Iterate over the planes until we reach the first plane (one loop)
    for (std::size_t i_p = 0u; i_p < n_planes; i_p++) {

        // Next surface index (circular)
        const std::size_t next_index = i_p == n_planes - 1 ? 0u : i_p + 1;

        // Helix on the current surface
        detail::helix<transform3> hlx(free_trk, &B);

        // bound vector on the current surface
        const auto bound_vec =
            c2.free_to_bound_vector(trfs[i_p], free_trk.vector());

        // Bound-to-free jacobian on the current surface
        const typename cartesian2<transform3>::bound_to_free_matrix
            bound_to_free_jacobi =
                c2.bound_to_free_jacobian(trfs[i_p], rectangle, bound_vec);

        // Get the intersection on the next surface
        const auto is = hpi(hlx, dindex_invalid, rectangle, trfs[next_index]);
        ASSERT_TRUE(is.status == intersection::status::e_inside);

        // Helical path length between two surfaces
        const auto path_length = is.path;

        // Add the path length to the total path length
        total_path_length += path_length;

        // Free transport jacobian between two surfaces
        const cartesian2<transform3>::free_matrix transport_jacobi =
            hlx.jacobian(path_length);

        // Reset the free track parameters with the position and direction on
        // the next surface
        free_trk.set_pos(is.p3);
        free_trk.set_dir(hlx.dir(path_length));

        // dtds
        const vector3 dtds = free_trk.qop() * vector::cross(free_trk.dir(), B);

        // Path correction
        const cartesian2<transform3>::free_matrix path_correction =
            c2.path_correction(free_trk.pos(), free_trk.dir(), dtds,
                               trfs[next_index], rectangle);

        // Correction term for the path variation
        const cartesian2<transform3>::free_matrix correction_term =
            matrix_operator().template identity<e_free_size, e_free_size>() +
            path_correction;

        // Free to bound jacobian
        const typename cartesian2<transform3>::free_to_bound_matrix
            free_to_bound_jacobi = c2.free_to_bound_jacobian(
                trfs[next_index], rectangle, free_trk.vector());

        // Full jacobian
        const typename cartesian2<transform3>::bound_matrix full_jacobi =
            free_to_bound_jacobi * correction_term * transport_jacobi *
            bound_to_free_jacobi;

        // Update the covariance at the next surface
        bound_cov =
            full_jacobi * bound_cov * matrix_operator().transpose(full_jacobi);
    }

    // Check if the total path length is the expected value
    ASSERT_NEAR(total_path_length, S, tolerance);

    // Check if the same covariance is obtained after one loop
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov, i, j), tolerance);
        }
    }
}
