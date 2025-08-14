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

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(propagator_state_t& propagation) const {

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

        // Set the noise to be expected by the navigator after all of the actors
        // are done and the covariance is up to date

        // Rough estimation of the track displacement at the next
        // surface
        // TODO: Find a better model
        constexpr scalar_type n_stddev{3.f};

        const auto& cov = bound_params.covariance();
        const scalar_type delta_phi{
            n_stddev *
            math::sqrt(getter::element(cov, e_bound_phi, e_bound_phi))};
        const scalar_type delta_theta{
            n_stddev *
            math::sqrt(getter::element(cov, e_bound_theta, e_bound_theta))};

        const scalar_type theta{bound_params.theta()};

        // Distance in spherical coordinates
        const scalar_type dist_2{
            math::sin(theta) * math::sin(theta + delta_theta) *
                math::cos(delta_phi) +
            math::cos(theta) * math::cos(theta + delta_theta)};

        navigation.set_external_tol(constant<scalar_type>::sqrt2 *
                                    math::sqrt(dist_2) * navigation());
    }
};

}  // namespace detray
