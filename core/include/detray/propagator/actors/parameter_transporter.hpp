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
struct parameter_transporter : actor {

    struct state {};

    /// Mask store visitor
    struct kernel {

        /// @name Type definitions for the struct
        /// @{

        // Transformation matching this struct
        using transform3_type = transform3_t;
        // scalar_type
        using scalar_type = typename transform3_type::scalar_type;
        // size type
        using size_type = typename transform3_type::size_type;
        // Matrix actor
        using matrix_operator = typename transform3_t::matrix_actor;
        // 2D matrix type
        template <size_type ROWS, size_type COLS>
        using matrix_type =
            typename matrix_operator::template matrix_type<ROWS, COLS>;
        // Shorthand vector/matrix types related to bound track parameters.
        using bound_vector = matrix_type<e_bound_size, 1>;
        using bound_matrix = matrix_type<e_bound_size, e_bound_size>;
        // Mapping from bound track parameters.
        using bound_to_free_matrix = matrix_type<e_free_size, e_bound_size>;
        // Shorthand vector/matrix types related to free track parameters.
        using free_vector = matrix_type<e_free_size, 1>;
        using free_matrix = matrix_type<e_free_size, e_free_size>;
        // Mapping from free track parameters.
        using free_to_bound_matrix = matrix_type<e_bound_size, e_free_size>;
        using free_to_path_matrix = matrix_type<1, e_free_size>;
        // Track helper
        using track_helper = detail::track_helper<matrix_operator>;
        // Vector in 3D space
        using vector3 = typename transform3_t::vector3;

        /// @}

        template <typename mask_group_t, typename index_t,
                  typename propagator_state_t>
        DETRAY_HOST_DEVICE inline void operator()(
            const mask_group_t& mask_group, const index_t& index,
            const transform3_type& trf3, propagator_state_t& propagation) {

            // Stepper and Navigator states
            auto& stepping = propagation._stepping;

            // Mask
            const auto& mask = mask_group[index];
            auto local_coordinate = mask.local_frame();

            // Free vector
            const auto& free_vec = stepping().vector();

            // Convert free to bound vector
            stepping._bound_params.set_vector(
                local_coordinate.free_to_bound_vector(trf3, free_vec));

            // Free to bound jacobian at the destination surface
            const free_to_bound_matrix free_to_bound_jacobian =
                local_coordinate.free_to_bound_jacobian(trf3, free_vec);

            // Transport jacobian in free coordinate
            free_matrix& free_transport_jacobian = stepping._jac_transport;

            // Path correction factor
            free_matrix path_correction = local_coordinate.path_correction(
                stepping().pos(), stepping().dir(), stepping.dtds(), trf3);

            const free_matrix correction_term =
                matrix_operator()
                    .template identity<e_free_size, e_free_size>() +
                path_correction;

            bound_matrix new_cov =
                matrix_operator().template zero<e_bound_size, e_bound_size>();

            if (propagation.param_type() == parameter_type::e_free) {

                const matrix_type<e_bound_size, e_free_size> full_jacobian =
                    free_to_bound_jacobian * correction_term *
                    free_transport_jacobian;

                new_cov = full_jacobian * stepping().covariance() *
                          matrix_operator().transpose(full_jacobian);

                propagation.set_param_type(parameter_type::e_bound);

            } else if (propagation.param_type() == parameter_type::e_bound) {
                // Bound to free jacobian at the departure surface
                const bound_to_free_matrix& bound_to_free_jacobian =
                    stepping._jac_to_global;

                stepping._full_jacobian =
                    free_to_bound_jacobian * correction_term *
                    free_transport_jacobian * bound_to_free_jacobian;

                new_cov = stepping._full_jacobian *
                          stepping._bound_params.covariance() *
                          matrix_operator().transpose(stepping._full_jacobian);
            }

            // Calculate surface-to-surface covariance transport
            stepping._bound_params.set_covariance(new_cov);
        }
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& /*actor_state*/,
                                       propagator_state_t& propagation) const {
        auto& navigation = propagation._navigation;

        // Do covariance transport when the track is on surface
        if (navigation.is_on_module()) {
            auto& stepping = propagation._stepping;

            const auto* det = navigation.detector();
            const auto& trf_store = det->transform_store();
            const auto& mask_store = det->mask_store();

            // Surface
            const auto& surface = navigation.current()->surface;

            mask_store.template visit<kernel>(
                surface.mask(), trf_store[surface.transform()], propagation);

            // Set surface link
            stepping._bound_params.set_surface_link(
                navigation.current_object());
        }
    }
};  // namespace detray

}  // namespace detray