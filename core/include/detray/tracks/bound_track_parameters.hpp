/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/geometry/barcode.hpp"
#include "detray/tracks/detail/track_helper.hpp"

namespace detray {

template <typename algebra_t>
struct bound_track_parameters {

    /// @name Type definitions for the struct
    /// @{
    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using point2_type = dpoint2D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using matrix_operator = dmatrix_operator<algebra_t>;

    // Shorthand vector/matrix types related to bound track parameters.
    using vector_type = bound_vector<algebra_t>;
    using covariance_type = bound_matrix<algebra_t>;

    // Track helper
    using track_helper = detail::track_helper<matrix_operator>;

    /// @}

    DETRAY_HOST_DEVICE
    bound_track_parameters()
        : m_barcode(),
          m_vector(matrix_operator().template zero<e_bound_size, 1>()),
          m_covariance(
              matrix_operator().template zero<e_bound_size, e_bound_size>()) {}

    DETRAY_HOST_DEVICE
    bound_track_parameters(const geometry::barcode sf_idx,
                           const vector_type& vec, const covariance_type& cov)
        : m_barcode(sf_idx), m_vector(vec), m_covariance(cov) {}

    /** @param rhs is the left hand side params for comparison
     **/
    DETRAY_HOST_DEVICE
    bool operator==(const bound_track_parameters& rhs) const {
        if (m_barcode != rhs.surface_link()) {
            return false;
        }

        for (unsigned int i = 0u; i < e_bound_size; i++) {
            const auto lhs_val = matrix_operator().element(m_vector, i, 0u);
            const auto rhs_val = matrix_operator().element(rhs.vector(), i, 0u);

            if (math::fabs(lhs_val - rhs_val) >
                std::numeric_limits<scalar_type>::epsilon()) {
                return false;
            }
        }
        for (unsigned int i = 0u; i < e_bound_size; i++) {
            for (unsigned int j = 0u; j < e_bound_size; j++) {
                const auto lhs_val =
                    matrix_operator().element(m_covariance, i, j);
                const auto rhs_val =
                    matrix_operator().element(rhs.covariance(), i, j);

                if (math::fabs(lhs_val - rhs_val) >
                    std::numeric_limits<scalar_type>::epsilon()) {
                    return false;
                }
            }
        }
        return true;
    }

    DETRAY_HOST_DEVICE
    const geometry::barcode& surface_link() const { return m_barcode; }

    DETRAY_HOST_DEVICE
    void set_surface_link(geometry::barcode link) { m_barcode = link; }

    DETRAY_HOST_DEVICE
    vector_type& vector() { return m_vector; }

    DETRAY_HOST_DEVICE
    const vector_type& vector() const { return m_vector; }

    DETRAY_HOST_DEVICE
    void set_vector(const vector_type& v) { m_vector = v; }

    DETRAY_HOST_DEVICE
    covariance_type& covariance() { return m_covariance; }

    DETRAY_HOST_DEVICE
    const covariance_type& covariance() const { return m_covariance; }

    DETRAY_HOST_DEVICE
    void set_covariance(const covariance_type& c) { m_covariance = c; }

    DETRAY_HOST_DEVICE
    point2_type bound_local() const {
        return track_helper().bound_local(m_vector);
    }

    DETRAY_HOST_DEVICE
    scalar_type phi() const {
        return matrix_operator().element(m_vector, e_bound_phi, 0u);
    }

    DETRAY_HOST_DEVICE
    scalar_type theta() const {
        return matrix_operator().element(m_vector, e_bound_theta, 0u);
    }

    DETRAY_HOST_DEVICE
    vector3_type dir() const { return track_helper().dir(m_vector); }

    DETRAY_HOST_DEVICE
    scalar_type time() const { return track_helper().time(m_vector); }

    DETRAY_HOST_DEVICE
    scalar_type qop() const { return track_helper().qop(m_vector); }

    DETRAY_HOST_DEVICE
    void set_qop(const scalar_type qop) {
        track_helper().set_qop(m_vector, qop);
    }

    DETRAY_HOST_DEVICE
    scalar_type qopT() const { return track_helper().qopT(m_vector); }

    DETRAY_HOST_DEVICE
    scalar_type qopz() const { return track_helper().qopz(m_vector); }

    DETRAY_HOST_DEVICE
    scalar_type p(const scalar_type q) const {
        return track_helper().p(m_vector, q);
    }

    DETRAY_HOST_DEVICE
    vector3_type mom(const scalar_type q) const {
        return track_helper().mom(m_vector, q);
    }

    DETRAY_HOST_DEVICE
    scalar_type pT(const scalar_type q) const {
        assert(this->qop() != 0.f);
        return math::fabs(q / this->qop() * getter::perp(this->dir()));
    }

    DETRAY_HOST_DEVICE
    scalar_type pz(const scalar_type q) const {
        assert(this->qop() != 0.f);
        return math::fabs(q / this->qop() * this->dir()[2]);
    }

    private:
    geometry::barcode m_barcode;
    vector_type m_vector;
    covariance_type m_covariance;
};

}  // namespace detray
