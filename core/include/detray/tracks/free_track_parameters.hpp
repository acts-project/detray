/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/definitions/units.hpp"
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

        matrix_operator().element(m_vector, e_free_pos0, 0u) = pos[0];
        matrix_operator().element(m_vector, e_free_pos1, 0u) = pos[1];
        matrix_operator().element(m_vector, e_free_pos2, 0u) = pos[2];
        matrix_operator().element(m_vector, e_free_time, 0u) = time;

        scalar_type p = getter::norm(mom);
        auto mom_norm = vector::normalize(mom);
        matrix_operator().element(m_vector, e_free_dir0, 0u) = mom_norm[0];
        matrix_operator().element(m_vector, e_free_dir1, 0u) = mom_norm[1];
        matrix_operator().element(m_vector, e_free_dir2, 0u) = mom_norm[2];
        matrix_operator().element(m_vector, e_free_qoverp, 0u) = q / p;
    }

    /** @param rhs is the left hand side params for comparison
     **/
    DETRAY_HOST_DEVICE
    bool operator==(const free_track_parameters& rhs) const {
        for (unsigned int i = 0u; i < e_free_size; i++) {
            const auto lhs_val = matrix_operator().element(m_vector, i, 0u);
            const auto rhs_val = matrix_operator().element(rhs.vector(), i, 0u);

            if (std::abs(lhs_val - rhs_val) >
                std::numeric_limits<scalar_type>::epsilon()) {
                return false;
            }
        }
        for (unsigned int i = 0u; i < e_free_size; i++) {
            for (unsigned int j = 0u; j < e_free_size; j++) {
                const auto lhs_val =
                    matrix_operator().element(m_covariance, i, j);
                const auto rhs_val =
                    matrix_operator().element(rhs.covariance(), i, j);

                if (std::abs(lhs_val - rhs_val) >
                    std::numeric_limits<scalar_type>::epsilon()) {
                    return false;
                }
            }
        }
        return true;
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
    scalar_type time() const { return track_helper().time(m_vector); }

    DETRAY_HOST_DEVICE
    scalar_type charge() const { return track_helper().charge(m_vector); }

    DETRAY_HOST_DEVICE
    scalar_type qop() const { return track_helper().qop(m_vector); }

    DETRAY_HOST_DEVICE
    void set_qop(const scalar_type qop) {
        matrix_operator().element(m_vector, e_free_qoverp, 0u) = qop;
    }

    DETRAY_HOST_DEVICE
    scalar_type qopT() const { return track_helper().qopT(m_vector); }

    DETRAY_HOST_DEVICE
    scalar_type p() const { return track_helper().p(m_vector); }

    DETRAY_HOST_DEVICE
    vector3 mom() const { return track_helper().mom(m_vector); }

    DETRAY_HOST_DEVICE
    scalar_type pT() const {
        return std::abs(1.f / this->qop() * getter::perp(this->dir()));
    }

    private:
    vector_type m_vector = matrix_operator().template zero<e_free_size, 1>();
    covariance_type m_covariance =
        matrix_operator().template zero<e_free_size, e_free_size>();
    scalar_type m_overstep_tolerance{-100 * unit<scalar_type>::um};
};

}  // namespace detray
