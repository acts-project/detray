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
        int n_stddev{1};
        /// Add adaptive mask tolerance to navigation
        bool estimate_scattering_noise{true};
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(const state& resetter_state,
                                       propagator_state_t& propagation) const {

        using scalar_type = dscalar<algebra_t>;

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

        // Track direction is known precisely, don't adjust tolerances
        if (!resetter_state.estimate_scattering_noise) {
            return;
        }

        // Set the noise to be expected by the navigator after all of the actors
        // are done and the covariance is up to date

        // Rough estimation of the track displacement at the next surface
        const auto& cov = bound_params.covariance();
        const scalar_type delta_phi{
            static_cast<scalar_type>(resetter_state.n_stddev) *
            math::sqrt(getter::element(cov, e_bound_phi, e_bound_phi))};
        const scalar_type delta_theta{
            static_cast<scalar_type>(resetter_state.n_stddev) *
            math::sqrt(getter::element(cov, e_bound_theta, e_bound_theta))};

        const scalar_type theta{bound_params.theta()};

        // Distance in spherical coordinates
        scalar_type dist_2{1.f -
                           (math::sin(theta) * math::sin(theta + delta_theta) *
                                math::cos(-delta_phi) +
                            math::cos(theta) * math::cos(theta + delta_theta))};

        // Sometimes the dist_2 value drops below zero due to floating point
        // errors, so normalize it to zero in these cases
        dist_2 = math::signbit(dist_2) ? 0.f : dist_2;

        // Guess the portal envelope distance if there is no next target
        const scalar_type path{
            navigation.is_exhausted()
                ? 5.f * unit<scalar_type>::mm
                : math::fabs(std::as_const(navigation).target().path)};

        // Noise component from the angle covariances
        scalar_type scattering_noise{constant<scalar_type>::sqrt2 * path *
                                     math::sqrt(dist_2)};

        // Parametrized noise component that scales with the path length
        const scalar_type q{stepping.particle_hypothesis().charge()};
        scalar_type accumulated_noise{1.f * unit<scalar_type>::GeV /
                                      stepping().p(q) *
                                      resetter_state.accumulated_error *
                                      math::fabs(stepping.path_length())};

        // Accumulate error until next measurement
        if (!sf.is_sensitive()) {
            scattering_noise += navigation.external_tol();
        }

        navigation.set_external_tol(accumulated_noise + scattering_noise);
        assert(std::isfinite(navigation.external_tol()));
    }
};

}  // namespace detray
