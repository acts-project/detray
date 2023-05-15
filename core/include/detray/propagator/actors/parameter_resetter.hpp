/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/tracks/detail/track_helper.hpp"

namespace detray {

template <typename transform3_t>
struct parameter_resetter : actor {

    using scalar_type = typename transform3_t::scalar_type;

    struct state {};

    /// Mask store visitor
    struct kernel {

        /// @name Type definitions for the struct
        /// @{

        // Matrix actor
        using matrix_operator = typename transform3_t::matrix_actor;

        /// @}

        template <typename mask_group_t, typename index_t,
                  typename stepper_state_t>
        DETRAY_HOST_DEVICE inline void operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3_t& trf3, stepper_state_t& stepping) const {

            // Note: How is it possible with "range"???
            const auto& mask = mask_group[index];

            auto local_coordinate = mask.local_frame();

            // Reset the free vector
            stepping().set_vector(local_coordinate.bound_to_free_vector(
                trf3, mask, stepping._bound_params.vector()));

            // Reset the path length
            stepping._s = 0;

            // Reset jacobian coordinate transformation at the current surface
            stepping._jac_to_global = local_coordinate.bound_to_free_jacobian(
                trf3, mask, stepping._bound_params.vector());

            // Reset jacobian transport to identity matrix
            matrix_operator().set_identity(stepping._jac_transport);
        }
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& /*resetter_state*/,
                                       propagator_state_t& propagation) const {

        auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // Do covariance transport when the track is on surface
        if (navigation.is_on_module()) {

            auto* det = navigation.detector();
            const auto& trf_store = det->transform_store();
            const auto& mask_store = det->mask_store();

            // Surface
            const auto& surface = navigation.current()->surface;

            mask_store.template visit<kernel>(
                surface.mask(), trf_store[surface.transform()], stepping);
        }
    }
};

}  // namespace detray
