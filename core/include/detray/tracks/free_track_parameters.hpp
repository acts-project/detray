/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/track_parametrization.hpp"
#include "detray/definitions/units.hpp"
#include "detray/tracks/detail/track_helper.hpp"

namespace detray {

template <typename algebra_t>
struct free_track_parameters {

    /// @name Type definitions for the struct
    /// @{
    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;
    using matrix_operator = dmatrix_operator<algebra_t>;

    // Shorthand vector/matrix types related to free track parameters.
    using vector_type = free_vector<algebra_t>;

    // Track helper
    using track_helper = detail::track_helper<matrix_operator>;

    /// @}

    DETRAY_HOST_DEVICE
    free_track_parameters()
        : m_vector(matrix_operator().template zero<e_free_size, 1>()){};

    DETRAY_HOST_DEVICE
    free_track_parameters(const vector_type& vec) : m_vector(vec) {}

    DETRAY_HOST_DEVICE
    free_track_parameters(const point3_type& pos, const scalar_type time,
                          const vector3_type& mom, const scalar_type q) {

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

            if (math::fabs(lhs_val - rhs_val) >
                std::numeric_limits<scalar_type>::epsilon()) {
                return false;
            }
        }

        return true;
    }

    DETRAY_HOST_DEVICE
    const vector_type& vector() const { return m_vector; }

    DETRAY_HOST_DEVICE
    void set_vector(const vector_type& v) { m_vector = v; }

    DETRAY_HOST_DEVICE
    point3_type pos() const { return track_helper().pos(m_vector); }

    DETRAY_HOST_DEVICE
    void set_pos(const vector3_type& pos) {
        track_helper().set_pos(m_vector, pos);
    }

    DETRAY_HOST_DEVICE
    vector3_type dir() const { return track_helper().dir(m_vector); }

    DETRAY_HOST_DEVICE
    void set_dir(const vector3_type& dir) {
        track_helper().set_dir(m_vector, dir);
    }

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
    vector_type m_vector = matrix_operator().template zero<e_free_size, 1>();
};

}  // namespace detray
