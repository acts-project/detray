/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"

namespace detray {

template <concepts::algebra algebra_t>
struct parameter_resetter : actor {

    struct state {
        /// Percentage of total track path to assume as accumulated error
        float accumulated_error{0.001f};
        /// Number of standard deviations to assume to model the scattering
        /// noise
        int n_stddev{2};
        /// Estimate mask tolerance for navigation to for compensate scattering
        bool estimate_scattering_noise{true};
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(const state& resetter_state,
                                       propagator_state_t& propagation) const {

        using scalar_t = dscalar<algebra_t>;

        auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // Do covariance transport when the track is on surface
        if (!(navigation.is_on_sensitive() ||
              navigation.encountered_sf_material())) {
            return;
        }

        // Update free params after bound params were changed by actors
        const auto sf = navigation.get_surface();
        const auto& bound_params = stepping.bound_params();
        stepping() =
            sf.bound_to_free_vector(propagation._context, bound_params);
        assert(!stepping().is_invalid());

        // Reset jacobian transport to identity matrix
        stepping.reset_transport_jacobian();

        // Track pos/dir is not known precisely: adjust navigation tolerances
        if (resetter_state.estimate_scattering_noise) {
            estimate_external_mask_tolerance(
                bound_params, propagation,
                static_cast<scalar_t>(resetter_state.n_stddev),
                resetter_state.accumulated_error);
        }
    }

    private:
    /// Estimate the mask tolerance that is needed for the navigation to include
    /// the next surface from the current covariance
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE constexpr void estimate_external_mask_tolerance(
        const bound_track_parameters<algebra_t>& bound_params,
        propagator_state_t& propagation, const dscalar<algebra_t> n_stddev,
        const dscalar<algebra_t> accumulated_error = 0.f) const {

        using scalar_t = dscalar<algebra_t>;

        auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // Set the noise to be expected by the navigator after all of the actors
        // are done and the covariance is up to date

        // Positional error on the current surface as an estimate of the error
        // on the next surface
        const auto& cov = bound_params.covariance();

        const scalar_t var_loc0{
            getter::element(cov, e_bound_loc0, e_bound_loc0)};
        const scalar_t var_loc1{
            getter::element(cov, e_bound_loc1, e_bound_loc1)};

        // Rough estimation of the track displacement at the next surface
        const scalar_t delta_phi{
            n_stddev *
            math::sqrt(getter::element(cov, e_bound_phi, e_bound_phi))};
        const scalar_t delta_theta{
            n_stddev *
            math::sqrt(getter::element(cov, e_bound_theta, e_bound_theta))};

        // Calculate the difference in cartesian coordinates, as the conversion
        // uses less trigonometric functions than calculating the distance in
        // spherical coordinates
        const scalar_t phi_err{bound_params.phi() + delta_phi};
        const scalar_t theta_err{bound_params.theta() + delta_theta};
        const scalar_t sin_theta_err{math::sin(theta_err)};

        dvector3D<algebra_t> displ{math::cos(phi_err) * sin_theta_err,
                                   math::sin(phi_err) * sin_theta_err,
                                   math::cos(theta_err)};

        // Guess the portal envelope distance if there is no next target
        const scalar_t path{
            navigation.is_exhausted()
                ? 5.f * unit<scalar_t>::mm
                : math::fabs(std::as_const(navigation).target().path())};

        displ = path * (displ - stepping().dir());

        // Parametrized noise component that scales with the path length
        // Accounts for material/interaction mismodelling
        const scalar_t q{stepping.particle_hypothesis().charge()};
        const scalar_t accumulated_noise{1.f * unit<scalar_t>::GeV /
                                         stepping().p(q) * accumulated_error *
                                         stepping.path_length()};

        navigation.set_external_tol(math::sqrt(
            n_stddev * n_stddev * (var_loc0 + var_loc1) +
            vector::dot(displ, displ) + accumulated_noise * accumulated_noise));

        // Clip to 5mm if the covariances are very large
        constexpr auto max_tol{5.f * unit<scalar_t>::mm};
        navigation.set_external_tol(navigation.external_tol() > max_tol
                                        ? max_tol
                                        : navigation.external_tol());

        assert(std::isfinite(navigation.external_tol()));
    }
};

}  // namespace detray
