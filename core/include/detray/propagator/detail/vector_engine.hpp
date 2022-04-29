/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parameterization.hpp"

namespace detray {

namespace detail {

template <typename scalar_t>
struct vector_engine {

    using vector3 = __plugin::vector3<scalar_t>;
    using point3 = __plugin::point3<scalar_t>;
    using size_type = __plugin::size_type;
    using transform3 = __plugin::transform3<scalar_t>;

    template <size_type ROWS, size_type COLS>
    using matrix_type = __plugin::matrix_type<scalar_t, ROWS, COLS>;
    using matrix_operator = standard_matrix_operator<scalar_t>;

    /** Function to get the position from free vector
     *
     * @param free_vec is the input free vector
     * @returns point3 position
     */
    DETRAY_HOST_DEVICE point3 pos(const free_vector& free_vec) const;

    /** Function to set the position of free vector
     *
     * @param free_vec is the input free vector
     */
    DETRAY_HOST_DEVICE void set_pos(free_vector& free_vec, const point3& pos);

    /** Function to get the direction from free vector
     *
     * @param free_vec is the input free vector
     * @returns vector3 direction
     */
    DETRAY_HOST_DEVICE vector3 dir(const free_vector& free_vec) const;

    /** Function to set the direction of free vector
     *
     * @param free_vec is the input free vector
     */
    DETRAY_HOST_DEVICE void set_dir(free_vector& free_vec, const vector3& dir);

    /** Function to get the local position of bound vector
     *
     * @param bound_vec is the input bound vector
     * @returns point3 local position
     */
    DETRAY_HOST_DEVICE point3 local(const bound_vector& bound_vec) const;

    /** Function to get the direction from bound vector
     *
     * @param bound_vec is the input bound vector
     * @returns vector3 direction
     */
    DETRAY_HOST_DEVICE vector3 dir(const bound_vector& bound_vec) const;

    /** Function to get the free vector from bound vector
     *
     * @param trf3 is the transform matrix
     * @param free_vec is the input free vector
     * @returns bound vector
     */
    DETRAY_HOST_DEVICE bound_vector free_to_bound_vector(
        const transform3& trf, const free_vector& free_vec) const;

    /** Function to get the bound vector from free vector
     *
     * @param trf3 is the transform matrix
     * @param bound_vec is the input bound vector
     * @returns free vector
     */
    DETRAY_HOST_DEVICE free_vector bound_to_free_vector(
        const transform3& trf, const bound_vector& bound_vec) const;
};

}  // namespace detail

}  // namespace detray

#include "detray/propagator/detail/vector_engine.ipp"