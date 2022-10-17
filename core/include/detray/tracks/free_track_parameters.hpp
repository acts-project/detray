/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/tracks/detail/track_helper.hpp"

namespace detray {

template <typename transform3_t>
struct free_track_parameters {

    /// @name Type definitions for the struct
    /// @{

    using transform3_type = transform3_t;
    using matrix_operator = typename transform3_type::matrix_actor;
    using size_type = typename transform3_type::size_type;
    using scalar_type = typename transform3_type::scalar_type;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename transform3_type::matrix_actor::template matrix_type<ROWS,
                                                                     COLS>;
    using vector3 = typename transform3_type::vector3;
    using point3 = typename transform3_type::point3;
    using point2 = typename transform3_type::point2;

    // Shorthand vector/matrix types related to free track parameters.
    using vector_type = matrix_type<e_free_size, 1>;
    using covariance_type = matrix_type<e_free_size, e_free_size>;

    // Track helper
    using track_helper = detail::track_helper<matrix_operator>;

    /// @}

    DETRAY_HOST_DEVICE
    free_track_parameters()
        : m_vector(matrix_operator().template zero<e_free_size, 1>()),
          m_covariance(
              matrix_operator().template zero<e_free_size, e_free_size>()){};

    DETRAY_HOST_DEVICE
    free_track_parameters(const vector_type& vec, const covariance_type& cov)
        : m_vector(vec), m_covariance(cov) {}

    DETRAY_HOST_DEVICE
    free_track_parameters(const point3& pos, const scalar_type time,
                          const vector3& mom, const scalar_type q) {

        matrix_operator().element(m_vector, e_free_pos0, 0) = pos[0];
        matrix_operator().element(m_vector, e_free_pos1, 0) = pos[1];
        matrix_operator().element(m_vector, e_free_pos2, 0) = pos[2];
        matrix_operator().element(m_vector, e_free_time, 0) = time;

        scalar_type p = getter::norm(mom);
        auto mom_norm = vector::normalize(mom);
        matrix_operator().element(m_vector, e_free_dir0, 0) = mom_norm[0];
        matrix_operator().element(m_vector, e_free_dir1, 0) = mom_norm[1];
        matrix_operator().element(m_vector, e_free_dir2, 0) = mom_norm[2];
        matrix_operator().element(m_vector, e_free_qoverp, 0) = q / p;
    }

    DETRAY_HOST_DEVICE
    const vector_type& vector() const { return m_vector; }

    DETRAY_HOST_DEVICE
    void set_vector(const vector_type& v) { m_vector = v; }

    DETRAY_HOST_DEVICE
    const covariance_type& covariance() const { return m_covariance; }

    DETRAY_HOST_DEVICE
    void set_covariance(const covariance_type& c) { m_covariance = c; }

    DETRAY_HOST_DEVICE
    scalar_type overstep_tolerance() const { return m_overstep_tolerance; }

    DETRAY_HOST_DEVICE
    void set_overstep_tolerance(const scalar_type tolerance) {
        m_overstep_tolerance = tolerance;
    }

    DETRAY_HOST_DEVICE
    point3 pos() const { return track_helper().pos(m_vector); }

    DETRAY_HOST_DEVICE
    void set_pos(const vector3& pos) { track_helper().set_pos(m_vector, pos); }

    DETRAY_HOST_DEVICE
    vector3 dir() const { return track_helper().dir(m_vector); }

    DETRAY_HOST_DEVICE
    void set_dir(const vector3& dir) { track_helper().set_dir(m_vector, dir); }

    DETRAY_HOST_DEVICE
    vector3 mom() const { return 1. / std::abs(this->qop()) * this->dir(); }

    DETRAY_HOST_DEVICE
    scalar_type time() const {
        return matrix_operator().element(m_vector, e_free_time, 0);
    }

    DETRAY_HOST_DEVICE
    scalar_type charge() const {
        return matrix_operator().element(m_vector, e_free_qoverp, 0) < 0 ? -1.
                                                                         : 1.;
    }

    DETRAY_HOST_DEVICE
    scalar_type qop() const {
        return matrix_operator().element(m_vector, e_free_qoverp, 0);
    }

    DETRAY_HOST_DEVICE
    scalar_type pT() const {
        auto dir = this->dir();
        return std::abs(1. / this->qop() * getter::perp(dir));
    }

    DETRAY_HOST_DEVICE
    void flip() {
        matrix_operator().element(m_vector, e_free_dir0, 0) *= -1.;
        matrix_operator().element(m_vector, e_free_dir1, 0) *= -1.;
        matrix_operator().element(m_vector, e_free_dir2, 0) *= -1.;
        matrix_operator().element(m_vector, e_free_qoverp, 0) *= -1.;
    }

    private:
    vector_type m_vector;
    covariance_type m_covariance;
    scalar_type m_overstep_tolerance = -1e-4;
};

}  // namespace detray