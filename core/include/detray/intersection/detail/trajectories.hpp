/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/matrix_helper.hpp"

// System include(s).
#include <climits>
#include <cmath>

namespace detray::detail {

/// @brief describes a straight-line trajectory
template <typename transform3_t>
class ray {
    public:
    using transform3_type = transform3_t;
    using scalar_type = typename transform3_type::scalar_type;
    using vector3 = typename transform3_type::vector3;
    using point3 = typename transform3_type::point3;

    /// Parametrized constructor that complies with track interface
    ///
    /// @param track the track state that should be approximated
    DETRAY_HOST_DEVICE ray(const free_track_parameters<transform3_type> &track)
        : _pos{track.pos()},
          _dir{track.dir()},
          _overstep_tolerance{track.overstep_tolerance()} {}

    /// Parametrized constructor that complies with track interface
    ///
    /// @param pos the track position
    /// @param dir the track momentum direction
    DETRAY_HOST_DEVICE
    ray(const point3 pos, const scalar_type /*time*/, const vector3 dir,
        const scalar_type /*q*/)
        : _pos{pos}, _dir{vector::normalize(dir)} {}

    /// @returns position on the ray (compatible with tracks/intersectors)
    DETRAY_HOST_DEVICE point3 pos() const { return _pos; }

    /// @returns position on the ray paramterized by path length
    DETRAY_HOST_DEVICE point3 pos(const scalar_type s) const {
        // Direction is always normalized in the constructor
        return _pos + s * _dir;
    }

    /// @param position new position on the ray
    DETRAY_HOST_DEVICE void set_pos(point3 pos) { _pos = pos; }

    /// @returns direction of the ray (compatible with tracks/intersectors)
    DETRAY_HOST_DEVICE vector3 dir() const { return _dir; }

    /// @returns direction of the ray paramterized by path length
    DETRAY_HOST_DEVICE vector3 dir(const scalar_type /*s*/) const {
        return this->dir();
    }

    /// @returns overstep tolerance to comply with track interface
    DETRAY_HOST_DEVICE
    scalar_type overstep_tolerance() const { return _overstep_tolerance; }

    /// Sets overstep tolerance to comply with track interface
    DETRAY_HOST_DEVICE
    void set_overstep_tolerance(const scalar_type tolerance) {
        _overstep_tolerance = tolerance;
    }

    private:
    /// origin of ray
    point3 _pos{0.f, 0.f, 0.f};
    /// direction of ray
    vector3 _dir{0.f, 0.f, 1.f};

    /// Overstep tolerance on a geometry surface
    scalar_type _overstep_tolerance{-1e-4f};
};

/// @brief describes a helical trajectory in a given B-field.
///
/// Helix class for the analytical solution of track propagation in
/// homogeneous B field. This Follows the notation of Eq (4.7) in
/// DOI:10.1007/978-3-030-65771-0
template <typename transform3_t>
class helix : public free_track_parameters<transform3_t> {
    public:
    using transform3_type = transform3_t;
    using scalar_type = typename transform3_type::scalar_type;
    using matrix_operator = typename transform3_type::matrix_actor;
    using vector3 = typename transform3_type::vector3;
    using point3 = typename transform3_type::point3;

    /// Free track parameters
    using free_track_parameters_type = free_track_parameters<transform3_t>;
    /// Size type
    using size_type = typename transform3_type::size_type;
    /// 2D Matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    using free_matrix = matrix_type<e_free_size, e_free_size>;
    using free_track_parameters_type::pos;
    using mat_helper = matrix_helper<matrix_operator>;

    DETRAY_HOST_DEVICE
    helix() = delete;

    /// Parametrized constructor
    ///
    /// @param pos the the origin of the helix
    /// @param time the time parameter
    /// @param dir the initial direction of momentum for the helix
    /// @param q the charge of the particle
    /// @param mag_field the magnetic field vector
    DETRAY_HOST_DEVICE
    helix(point3 pos, scalar_type time, vector3 dir, scalar_type q,
          vector3 const *const mag_field)
        : free_track_parameters_type(pos, time, dir, q),
          _mag_field(mag_field) {}

    /// Parametrized constructor
    ///
    /// @param vertex the underlying track parametrization
    /// @param mag_fied the magnetic field vector
    DETRAY_HOST_DEVICE
    helix(const free_track_parameters_type track,
          vector3 const *const mag_field)
        : free_track_parameters_type(track), _mag_field(mag_field) {

        // Normalized B field
        _h0 = vector::normalize(*_mag_field);

        // Normalized tangent vector
        _t0 = vector::normalize(free_track_parameters_type::mom());

        // Normalized _h0 X _t0
        _n0 = vector::normalize(vector::cross(_h0, _t0));

        // Magnitude of _h0 X _t0
        _alpha = getter::norm(vector::cross(_h0, _t0));

        // Dot product of _h0 X _t0
        _delta = vector::dot(_h0, _t0);

        // Path length scaler
        _K = -1.f * free_track_parameters_type::qop() *
             getter::norm(*_mag_field);

        // Get longitudinal momentum parallel to B field
        scalar_type pz = vector::dot(free_track_parameters_type::mom(), _h0);

        // Get transverse momentum perpendicular to B field
        vector3 pT = free_track_parameters_type::mom() - pz * _h0;

        // R [mm] =  pT [GeV] / B [T] in natrual unit
        _R = getter::norm(pT) / getter::norm(*_mag_field);

        // Handle the case of pT ~ 0
        if (getter::norm(pT) < 1e-6f) {
            _vz_over_vt = std::numeric_limits<scalar_type>::infinity();
        } else {
            // Get vz over vt in new coordinate
            _vz_over_vt = pz / getter::norm(pT);
        }
    }

    /// @returns the radius of helix
    DETRAY_HOST_DEVICE
    scalar_type radius() const { return _R; }

    /// @returns the position after propagating the path length of s
    DETRAY_HOST_DEVICE
    point3 operator()(const scalar_type s) const { return this->pos(s); }

    /// @returns the position after propagating the path length of s
    DETRAY_HOST_DEVICE
    point3 pos(const scalar_type s) const {

        // Handle the case of pT ~ 0
        if (_vz_over_vt == std::numeric_limits<scalar_type>::infinity()) {
            return free_track_parameters_type::pos() + s * _h0;
        }

        point3 ret = free_track_parameters_type::pos();
        ret = ret + _delta / _K * (_K * s - std::sin(_K * s)) * _h0;
        ret = ret + std::sin(_K * s) / _K * _t0;
        ret = ret + _alpha / _K * (1.f - std::cos(_K * s)) * _n0;

        return ret;
    }

    /// @returns the tangential vector after propagating the path length of s
    DETRAY_HOST_DEVICE
    vector3 dir(const scalar_type s) const {

        // Handle the case of pT ~ 0
        if (_vz_over_vt == std::numeric_limits<scalar_type>::infinity()) {
            return free_track_parameters_type::dir();
        }

        vector3 ret{0.f, 0.f, 0.f};

        ret = ret + _delta * (1 - std::cos(_K * s)) * _h0;
        ret = ret + std::cos(_K * s) * _t0;
        ret = ret + _alpha * std::sin(_K * s) * _n0;

        return ret;
    }

    /// @returns the transport jacobian after propagating the path length of s
    DETRAY_HOST_DEVICE
    free_matrix jacobian(const scalar_type s) const {

        free_matrix ret =
            matrix_operator().template zero<e_free_size, e_free_size>();

        const matrix_type<3, 3> I33 =
            matrix_operator().template identity<3, 3>();
        const matrix_type<3, 3> Z33 = matrix_operator().template zero<3, 3>();

        // Notations
        // r = position
        // t = direction
        // l = qoverp

        // Get drdr
        auto drdr = I33;
        matrix_operator().set_block(ret, drdr, e_free_pos0, e_free_pos0);

        // Get dtdr
        auto dtdr = Z33;
        matrix_operator().set_block(ret, dtdr, e_free_dir0, e_free_pos0);

        // Get drdt
        auto drdt = Z33;

        drdt = drdt + std::sin(_K * s) / _K * I33;

        const auto H0 = mat_helper().column_wise_multiply(I33, _h0);
        drdt = drdt + (_K * s - std::sin(_K * s)) / _K *
                          mat_helper().column_wise_multiply(
                              matrix_operator().transpose(H0), _h0);

        drdt = drdt + (std::cos(_K * s) - 1.f) / _K *
                          mat_helper().column_wise_cross(I33, _h0);

        matrix_operator().set_block(ret, drdt, e_free_pos0, e_free_dir0);

        // Get dtdt
        auto dtdt = Z33;
        dtdt = dtdt + std::cos(_K * s) * I33;
        dtdt = dtdt + (1.f - std::cos(_K * s)) *
                          mat_helper().column_wise_multiply(
                              matrix_operator().transpose(H0), _h0);
        dtdt =
            dtdt - std::sin(_K * s) * mat_helper().column_wise_cross(I33, _h0);

        matrix_operator().set_block(ret, dtdt, e_free_dir0, e_free_dir0);

        // Get drdl
        vector3 drdl = 1.f / free_track_parameters_type::qop() *
                       (s * this->dir(s) + free_track_parameters_type::pos() -
                        this->pos(s));

        matrix_operator().set_block(ret, drdl, e_free_pos0, e_free_qoverp);

        // Get dtdl
        vector3 dtdl =
            _alpha * _K * s / free_track_parameters_type::qop() * _n0;

        matrix_operator().set_block(ret, dtdl, e_free_dir0, e_free_qoverp);

        // 3x3 and 7x7 element is 1 (Maybe?)
        matrix_operator().element(ret, e_free_time, e_free_time) = 1.f;
        matrix_operator().element(ret, e_free_qoverp, e_free_qoverp) = 1.f;

        return ret;
    }

    /// B field
    vector3 const *_mag_field;

    /// Normalized b field
    vector3 _h0;

    /// Normalized tangent vector
    vector3 _t0;

    /// Normalized _h0 X _t0
    vector3 _n0;

    /// Magnitude of _h0 X _t0
    scalar_type _alpha;

    /// Dot product of _h0 X _t0
    scalar_type _delta;

    /// Path length scaler
    scalar_type _K;

    /// Radius [mm] of helix
    scalar_type _R;

    /// Velocity in new z axis divided by transverse velocity
    scalar_type _vz_over_vt;
};

}  // namespace detray::detail
