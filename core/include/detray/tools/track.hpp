/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// detray tools
#include "detray/definitions/qualifiers.hpp"
#include "detray/tools/track_parameterization.hpp"
#include "detray/utils/indexing.hpp"

namespace detray {

struct bound_track_parameters {

    using point2 = __plugin::point2<scalar>;
    using point3 = __plugin::point3<scalar>;
    using vector2 = __plugin::vector2<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using vector_type = bound_vector;
    using covariance_type = bound_matrix;
    using jacobian_type = bound_matrix;

    bound_track_parameters(const dindex& sf_idx, const vector_type& params,
                           const covariance_type& cov)
        : _surface_link(sf_idx), _vector(params), _covariance(cov) {}

    DETRAY_HOST_DEVICE
    vector_type& vector() { return _vector; }

    DETRAY_HOST_DEVICE
    covariance_type& covariance() { return _covariance; }

    DETRAY_HOST_DEVICE
    point2 local() const {
        return point2{getter::element(_vector, e_bound_loc0, 0),
                      getter::element(_vector, e_bound_loc1, 0)};
    }

    DETRAY_HOST_DEVICE
    scalar phi() const { return getter::element(_vector, e_bound_phi, 0); }

    DETRAY_HOST_DEVICE
    scalar theta() const { return getter::element(_vector, e_bound_theta, 0); }

    template <typename surface_container_t, typename transform_container_t>
    DETRAY_HOST_DEVICE point3
    pos(const typename transform_container_t::context& ctx,
        const surface_container_t& surfaces,
        const transform_container_t& trfs) const {

        auto tidx = surfaces[_surface_link].transform();
        const auto& trf = trfs.contextual_transform(ctx, tidx);
        point3 global = trf.point_to_global(
            point3{getter::element(_vector, e_bound_loc0, 0),
                   getter::element(_vector, e_bound_loc1, 0), 0.});
        return global;
    }

    template <typename detector_t>
    DETRAY_HOST_DEVICE vector3 pos(const typename detector_t::context& ctx,
                                   const detector_t& det) const {
        const auto& surfaces = det.surfaces();
        const auto& trfs = det.transforms(ctx);
        return this->pos(ctx, surfaces, trfs);
    }

    DETRAY_HOST_DEVICE
    vector3 dir() const {
        const auto& phi = getter::element(_vector, e_bound_phi, 0);
        const auto& theta = getter::element(_vector, e_bound_theta, 0);
        const auto cosTheta = std::cos(theta);
        const auto sinTheta = std::sin(theta);
        return {
            std::cos(phi) * sinTheta,
            std::sin(phi) * sinTheta,
            cosTheta,
        };
    }

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
    using point3 = __plugin::vector3<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using vector_type = free_vector;
    using covariance_type = free_sym_matrix;
    using jacobian_type = free_matrix;

    free_track_parameters(const vector_type& params, const covariance_type& cov)
        : _vector(params), _covariance(cov) {}

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
    vector_type& vector() { return _vector; }

    DETRAY_HOST_DEVICE
    covariance_type& covariance() { return _covariance; }

    DETRAY_HOST_DEVICE
    scalar overstep_tolerance() const { return _overstep_tolerance; }

    DETRAY_HOST_DEVICE
    void set_overstep_tolerance(scalar tolerance) {
        _overstep_tolerance = tolerance;
    }

    DETRAY_HOST_DEVICE
    point3 pos() const {
        return point3{getter::element(_vector, e_free_pos0, 0),
                      getter::element(_vector, e_free_pos1, 0),
                      getter::element(_vector, e_free_pos2, 0)};
    }

    DETRAY_HOST_DEVICE
    void set_pos(const vector3& pos) {
        getter::element(_vector, e_free_pos0, 0) = pos[0];
        getter::element(_vector, e_free_pos1, 0) = pos[1];
        getter::element(_vector, e_free_pos2, 0) = pos[2];
    }

    DETRAY_HOST_DEVICE
    vector3 dir() const {
        return vector3{getter::element(_vector, e_free_dir0, 0),
                       getter::element(_vector, e_free_dir1, 0),
                       getter::element(_vector, e_free_dir2, 0)};
    }

    DETRAY_HOST_DEVICE
    scalar time() const { return getter::element(_vector, e_free_time, 0); }

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

    private:
    vector_type _vector;
    covariance_type _covariance;
    scalar _overstep_tolerance = -1e-4;
};

}  // namespace detray
