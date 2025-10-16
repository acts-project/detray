/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
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
#include "detray/propagator/detail/noise_estimation.hpp"
#include "detray/propagator/propagation_config.hpp"
#include "detray/utils/log.hpp"

namespace detray {

template <concepts::algebra algebra_t>
struct parameter_resetter : actor {

    struct state {

        /// Default construction
        state() = default;

        /// Build from propagation configuration set by the user
        DETRAY_HOST_DEVICE
        explicit constexpr state(const propagation::config& cfg)
            : accumulated_error{cfg.navigation.accumulated_error},
              n_stddev{cfg.navigation.n_scattering_stddev},
              estimate_scattering_noise{
                  cfg.navigation.estimate_scattering_noise} {}

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
        DETRAY_VERBOSE_HOST_DEVICE("Update the free track parameters");

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
            detail::estimate_external_mask_tolerance(
                bound_params, propagation,
                static_cast<scalar_t>(resetter_state.n_stddev),
                resetter_state.accumulated_error);
        }
    }
};

}  // namespace detray
