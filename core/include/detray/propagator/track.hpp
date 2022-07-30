/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// detray tools
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/track_parameterization.hpp"
#include "detray/propagator/detail/vector_engine.hpp"

namespace detray {

struct bound_track_parameters {

    using point2 = __plugin::point2<scalar>;
    using point3 = __plugin::point3<scalar>;
    using vector2 = __plugin::vector2<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using vector_type = bound_vector;
    using covariance_type = bound_matrix;
    using jacobian_type = bound_matrix;
    using vector_engine = detail::vector_engine<scalar>;
    using matrix_operator = standard_matrix_operator<scalar>;

    DETRAY_HOST_DEVICE
    bound_track_parameters()
        : _vector(matrix_operator().template zero<e_bound_size, 1>()),
          _covariance(
              matrix_operator().template zero<e_bound_size, e_bound_size>()) {}

    DETRAY_HOST_DEVICE
    bound_track_parameters(const dindex& sf_idx, const vector_type& params,
                           const covariance_type& cov)
        : _surface_link(sf_idx), _vector(params), _covariance(cov) {}

    DETRAY_HOST_DEVICE
    const dindex& surface_link() const { return _surface_link; }

    DETRAY_HOST_DEVICE
    const vector_type& vector() const { return _vector; }

    DETRAY_HOST_DEVICE
    void set_vector(const vector_type& v) { _vector = v; }

    DETRAY_HOST_DEVICE
    const covariance_type& covariance() const { return _covariance; }

    DETRAY_HOST_DEVICE
    void set_covariance(const covariance_type& c) { _covariance = c; }

    DETRAY_HOST_DEVICE
    point3 local() const { return vector_engine().local(_vector); }

    DETRAY_HOST_DEVICE
    scalar phi() const { return getter::element(_vector, e_bound_phi, 0); }

    DETRAY_HOST_DEVICE
    scalar theta() const { return getter::element(_vector, e_bound_theta, 0); }

    DETRAY_HOST_DEVICE
    vector3 dir() const { return vector_engine().dir(_vector); }

    DETRAY_HOST_DEVICE
    scalar time() const { return getter::element(_vector, e_bound_time, 0); }

    DETRAY_HOST_DEVICE
    scalar charge() const {
        if (getter::element(_vector, e_bound_qoverp, 0) < 0) {
            return -1.;
        } else {
            return 1.;
        }
    }

    DETRAY_HOST_DEVICE
    scalar qop() const { return getter::element(_vector, e_bound_qoverp, 0); }

    private:
    dindex _surface_link;
    vector_type _vector;
    covariance_type _covariance;
};

struct free_track_parameters {
    using point2 = __plugin::vector2<scalar>;
    using point3 = __plugin::vector3<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using vector2 = __plugin::vector2<scalar>;
    using vector_type = free_vector;
    using covariance_type = free_sym_matrix;
    using jacobian_type = free_matrix;
    using vector_engine = detail::vector_engine<scalar>;
    using matrix_operator = standard_matrix_operator<scalar>;

    DETRAY_HOST_DEVICE
    free_track_parameters()
        : _vector(matrix_operator().template zero<e_free_size, 1>()),
          _covariance(
              matrix_operator().template zero<e_free_size, e_free_size>()){};

    DETRAY_HOST_DEVICE
    free_track_parameters(const vector_type& params, const covariance_type& cov)
        : _vector(params), _covariance(cov) {}

    DETRAY_HOST_DEVICE
    free_track_parameters(const point3& pos, const scalar& time,
                          const vector3& mom, const scalar& q) {

        getter::element(_vector, e_free_pos0, 0) = pos[0];
        getter::element(_vector, e_free_pos1, 0) = pos[1];
        getter::element(_vector, e_free_pos2, 0) = pos[2];
        getter::element(_vector, e_free_time, 0) = time;

        scalar p = getter::norm(mom);
        auto mom_norm = vector::normalize(mom);
        getter::element(_vector, e_free_dir0, 0) = mom_norm[0];
        getter::element(_vector, e_free_dir1, 0) = mom_norm[1];
        getter::element(_vector, e_free_dir2, 0) = mom_norm[2];
        getter::element(_vector, e_free_qoverp, 0) = q / p;
    }

    DETRAY_HOST_DEVICE
    const vector_type& vector() const { return _vector; }

    DETRAY_HOST_DEVICE
    void set_vector(const vector_type& v) { _vector = v; }

    DETRAY_HOST_DEVICE
    const covariance_type& covariance() const { return _covariance; }

    DETRAY_HOST_DEVICE
    void set_covariance(const covariance_type& c) { _covariance = c; }

    DETRAY_HOST_DEVICE
    scalar overstep_tolerance() const { return _overstep_tolerance; }

    DETRAY_HOST_DEVICE
    void set_overstep_tolerance(scalar tolerance) {
        _overstep_tolerance = tolerance;
    }

    DETRAY_HOST_DEVICE
    point3 pos() const { return vector_engine().pos(_vector); }

    DETRAY_HOST_DEVICE
    void set_pos(const vector3& pos) { vector_engine().set_pos(_vector, pos); }

    DETRAY_HOST_DEVICE
    vector3 dir() const { return vector_engine().dir(_vector); }

    DETRAY_HOST_DEVICE
    void set_dir(const vector3& dir) { vector_engine().set_dir(_vector, dir); }

    DETRAY_HOST_DEVICE
    vector3 mom() const { return 1. / std::abs(this->qop()) * this->dir(); }

    DETRAY_HOST_DEVICE
    scalar time() const { return getter::element(_vector, e_free_time, 0); }

    DETRAY_HOST_DEVICE
    void set_time(const scalar time) {
        getter::element(_vector, e_free_time, 0) = time;
    }

    DETRAY_HOST_DEVICE
    scalar charge() const {
        if (getter::element(_vector, e_free_qoverp, 0) < 0) {
            return -1.;
        } else {
            return 1.;
        }
    }

    DETRAY_HOST_DEVICE
    scalar qop() const { return getter::element(_vector, e_free_qoverp, 0); }

    DETRAY_HOST_DEVICE
    void set_qop(const scalar qop) {
        getter::element(_vector, e_free_qoverp, 0) = qop;
    }

    DETRAY_HOST_DEVICE
    scalar pT() const {
        auto dir = this->dir();
        return std::abs(1. / this->qop() *
                        getter::perp(vector2{dir[0], dir[1]}));
    }

    DETRAY_HOST_DEVICE
    void flip() {
        getter::element(_vector, e_free_dir0, 0) *= -1.;
        getter::element(_vector, e_free_dir1, 0) *= -1.;
        getter::element(_vector, e_free_dir2, 0) *= -1.;
        getter::element(_vector, e_free_qoverp, 0) *= -1.;
    }

    private:
    vector_type _vector;
    covariance_type _covariance;
    scalar _overstep_tolerance = -1e-4;
};

}  // namespace detray
