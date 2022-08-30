/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/definitions/qualifiers.hpp"

namespace detray {

namespace detail {

template <typename transform3_t>
struct bound_to_bound_covariance_update {

    /// @name Type definitions for the struct
    /// @{

    // Transformation matching this struct
    using transform3_type = transform3_t;
    // scalar_type
    using scalar_type = typename transform3_type::scalar_type;
    // size type
    using size_type = typename transform3_type::size_type;
    // Matrix actor
    using matrix_actor = typename transform3_t::matrix_actor;
    // 2D matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type = typename matrix_actor::template matrix_type<ROWS, COLS>;
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
    using track_helper = detail::track_helper<matrix_actor>;

    /// @}

    using output_type = bool;

    template <typename mask_group_t, typename propagator_state_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const mask_group_t& mask_group, propagator_state_t& propagation) const {

        // Stepper and Navigator states
        auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // Retrieve surfaces and transform store
        const auto& det = navigation.detector();
        const auto& surface_container = det->surfaces();
        const auto& tf_store = det->transform_store();

        // Intersection
        const auto& is = navigation.current();

        // Surface
        const auto& sf = surface_container[is->index];

        // Transform
        const auto& trf3 = tf_store[sf.transform()];

        // Mask
        const auto& mask = mask_group[is->mask_index];
        auto local_coordinate = mask.local_type();

        // Free vector
        const auto& free_vec = stepping().vector();

        // Convert free to bound vector
        const auto bound_vec =
            local_coordinate.free_to_bound_vector(trf3, free_vec);

        // Free to bound jacobian at the destination surface
        const free_to_bound_matrix free_to_bound_jacobian =
            local_coordinate.free_to_bound_jacobian(trf3, mask, free_vec);

        // Transport jacobian in free coordinate
        free_matrix& free_transport_jacobian = stepping._jac_transport;

        // Path correction factor

        // Bound to free jacobian at the departure surface
        const bound_to_free_matrix& bound_to_free_jacobian =
            stepping._jac_to_global;

        /*
        // Run over the masks belonged to the surface
        for (const auto& mask : range(mask_group, mask_range)) {

            auto local_coordinate = mask.local_type();

            local_coordinate.
        }
        */
        return true;
    }

    /** Function to get the full jacobian for surface-to-surface propagation
     *
     * @param trf3 is the transform matrix of the destination surface
     * @param mask is the mask of the destination surface
     * @param free_vec is the input free vector at the destination surface
     * @param bound_to_free_jacobian is the coordinate transform jacobian
     * at departure surface.
     * @param free_transport_jacobian is the transport jacobian
     * @param free_to_path_derivative is the derivative of free parameter w.r.t
     * path
     * @returns bound to bound jacobian
     */
    /*
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline bound_matrix bound_to_bound_jacobian(
        const transform3_t& trf3, const mask_t& mask,
        const free_vector& free_vec,
        const bound_to_free_matrix& bound_to_free_jacobian,
        const free_matrix& free_transport_jacobian,
        const free_matrix& path_correction) {

        using local_type = typename mask_t::local_type;

        const auto free_to_bound_matrix free_to_bound_jacobian =
            local_type().free_to_bound_jacobian(trf3, mask, free_vec);

        return free_to_bound_jacobian * free_transport_jacobian *
               bound_to_free_jacobian;
    }
    */
};

}  // namespace detail

}  // namespace detray