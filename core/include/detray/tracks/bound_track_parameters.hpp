/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/tracks/detail/track_helper.hpp"

namespace detray {

template <typename transform3_t>
struct bound_track_parameters {

    /// @name Type definitions for the struct
    /// @{

    using transform3_type = transform3_t;
    using matrix_actor = typename transform3_type::matrix_actor;
    using size_type = typename transform3_type::size_type;
    using scalar_type = typename transform3_type::scalar_type;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename transform3_type::matrix_actor::template matrix_type<ROWS,
                                                                     COLS>;
    using vector3 = typename transform3_type::vector3;
    using point3 = typename transform3_type::point3;
    using point2 = typename transform3_type::point2;

    // Shorthand vector/matrix types related to bound track parameters.
    using vector_type = matrix_type<e_bound_size, 1>;
    using covariance_type = matrix_type<e_bound_size, e_bound_size>;

    // Track helper
    using track_helper = detail::track_helper<matrix_actor>;

    /// @}

    DETRAY_HOST_DEVICE
    bound_track_parameters()
        : m_surface_link(dindex_invalid),
          m_vector(matrix_actor().template zero<e_bound_size, 1>()),
          m_covariance(
              matrix_actor().template zero<e_bound_size, e_bound_size>()) {}

    DETRAY_HOST_DEVICE
    bound_track_parameters(const std::size_t sf_idx, const vector_type& vec,
                           const covariance_type& cov)
        : m_surface_link(sf_idx), m_vector(vec), m_covariance(cov) {}

    DETRAY_HOST_DEVICE
    const std::size_t& surface_link() const { return m_surface_link; }

    DETRAY_HOST_DEVICE
    vector_type& vector() { return m_vector; }

    DETRAY_HOST_DEVICE
    const vector_type& vector() const { return m_vector; }

    DETRAY_HOST_DEVICE
    void set_vector(const vector_type& v) { m_vector = v; }

    DETRAY_HOST_DEVICE
    const covariance_type& covariance() const { return m_covariance; }

    DETRAY_HOST_DEVICE
    void set_covariance(const covariance_type& c) { m_covariance = c; }

    DETRAY_HOST_DEVICE
    point2 local() const { return track_helper().local(m_vector); }

    DETRAY_HOST_DEVICE
    scalar_type phi() const {
        return matrix_actor().element(m_vector, e_bound_phi, 0);
    }

    DETRAY_HOST_DEVICE
    scalar_type theta() const {
        return matrix_actor().element(m_vector, e_bound_theta, 0);
    }

    DETRAY_HOST_DEVICE
    vector3 dir() const { return track_helper().dir(m_vector); }

    DETRAY_HOST_DEVICE
    scalar_type time() const {
        return matrix_actor().element(m_vector, e_bound_time, 0);
    }

    DETRAY_HOST_DEVICE
    scalar_type charge() const {
        if (matrix_actor().element(m_vector, e_bound_qoverp, 0) < 0) {
            return -1.;
        } else {
            return 1.;
        }
    }

    DETRAY_HOST_DEVICE
    scalar_type qop() const {
        return matrix_actor().element(m_vector, e_bound_qoverp, 0);
    }

    private:
    std::size_t m_surface_link;
    vector_type m_vector;
    covariance_type m_covariance;
};

}  // namespace detray