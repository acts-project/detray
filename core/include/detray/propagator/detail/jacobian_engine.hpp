/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parameterization.hpp"
#include "detray/propagator/detail/vector_engine.hpp"

namespace detray {

namespace detail {

// Jacobian engine assuming "Planar Surface with Cartesian coordinate"
// NOTE: We may inherit jacobian engine to address different
// types of surface, if required
template <typename scalar_t>
struct jacobian_engine {

    using vector3 = __plugin::vector3<scalar_t>;
    template <__plugin::size_type ROWS, __plugin::size_type COLS>
    using matrix_type = __plugin::matrix_type<scalar_t, ROWS, COLS>;
    using vector_engine = detail::vector_engine<scalar_t>;
    using matrix_operator = typename vector_engine::matrix_operator;
    using transform3 = typename vector_engine::transform3;

    /** Function to get the jacobian for bound to free coordinate transform
     *
     * @param trf3 is the transform matrix
     * @param bound_vec is the input bound vector
     * @returns bound to free jacobian
     */
    DETRAY_HOST_DEVICE inline bound_to_free_matrix bound_to_free_coordinate(
        const transform3& trf3, const bound_vector& bound_vec) const;

    /** Function to get the jacobian for free to bound coordinate transform
     *
     * @param trf3 is the transform matrix
     * @param free_vec is the input free vector
     * @returns free to bound jacobian
     */
    DETRAY_HOST_DEVICE inline free_to_bound_matrix free_to_bound_coordinate(
        const transform3& trf3, const free_vector& free_vec) const;

    /** Function to calculate the path correction term for transport jacobian,
     * which is caused by the geometry constraint.
     *
     * @param trf3 is the transform matrix
     * @param free_vec is the input free vector
     * @returns free to path matrix
     */
    DETRAY_HOST_DEVICE inline free_to_path_matrix free_to_path_correction(
        const transform3& trf3, const free_vector& free_vec) const;
};

}  // namespace detail

}  // namespace detray

#include "detray/propagator/detail/jacobian_engine.ipp"