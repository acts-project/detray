/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// system include
#include <climits>
#include <cmath>

// detray include(s)
#include "detray/definitions/units.hpp"
#include "detray/propagator/track.hpp"
#include "detray/utils/algebra_helpers.hpp"

namespace detray {

namespace detail {

/// @brief describes a straight-line trajectory
class ray {
    public:
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /// Parametrized constructor that complies with track interface
    ///
    /// @param track the track state that should be approximated
    template <typename track_t>
    DETRAY_HOST_DEVICE ray(const track_t &track)
        : _pos(track.pos()), _dir(track.dir()) {}

    /// Parametrized constructor that complies with track interface
    ///
    /// @param pos the track position
    /// @param dir the track momentum direction
    DETRAY_HOST_DEVICE
    ray(const point3 pos, const scalar /*time*/, const vector3 dir,
        const scalar /*q*/)
        : _pos(pos), _dir(vector::normalize(dir)) {}

    /// @returns position on the ray (compatible with tracks/intersectors)
    DETRAY_HOST_DEVICE point3 pos() const { return _pos; }

    /// @returns position on the ray paramterized by path length
    DETRAY_HOST_DEVICE point3 pos(const scalar s) const {
        // Direction is always normalized in the constructor
        return _pos + s * _dir;
    }

    /// @param position new position on the ray
    DETRAY_HOST_DEVICE void set_pos(point3 pos) { _pos = pos; }

    /// @returns direction of the ray (compatible with tracks/intersectors)
    DETRAY_HOST_DEVICE vector3 dir() const { return _dir; }

    /// @returns direction of the ray paramterized by path length
    DETRAY_HOST_DEVICE vector3 dir(const scalar /*s*/) const {
        return this->dir();
    }

    /// @returns overstep tolerance to comply with track interface
    DETRAY_HOST_DEVICE
    scalar overstep_tolerance() const { return _overstep_tolerance; }

    private:
    /// origin of ray
    point3 _pos{0., 0., 0.};
    /// direction of ray
    vector3 _dir{0., 0., 1.};

    /// Overstep tolerance on a geometry surface
    scalar _overstep_tolerance{-1e-4};
};

/// @brief describes a helical trajectory in a given B-field.
///
/// Helix class for the analytical solution of track propagation in
/// homogeneous B field. This Follows the notation of Eq (4.7) in
/// DOI:10.1007/978-3-030-65771-0
class helix : public free_track_parameters {
    public:
    using point3 = __plugin::point3<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using vector2 = __plugin::vector2<scalar>;
    using matrix_operator = standard_matrix_operator<scalar>;
    using size_type = __plugin::size_type;
    template <size_type ROWS, size_type COLS>
    using matrix_type = __plugin::matrix_type<scalar, ROWS, COLS>;

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
    helix(point3 pos, scalar time, vector3 dir, scalar q,
          vector3 const *const mag_field)
        : free_track_parameters(pos, time, dir, q), _mag_field(mag_field) {}

    /// Parametrized constructor
    ///
    /// @param vertex the underlying track parametrization
    /// @param mag_fied the magnetic field vector
    DETRAY_HOST_DEVICE
    helix(const free_track_parameters vertex, vector3 const *const mag_field)
        : free_track_parameters(vertex), _mag_field(mag_field) {

        // Normalized B field
        _h0 = vector::normalize(*_mag_field);

        // Normalized tangent vector
        _t0 = vector::normalize(free_track_parameters::mom());

        // Normalized _h0 X _t0
        _n0 = vector::normalize(vector::cross(_h0, _t0));

        // Magnitude of _h0 X _t0
        _alpha = getter::norm(vector::cross(_h0, _t0));

        // Dot product of _h0 X _t0
        _delta = vector::dot(_h0, _t0);

        // Path length scaler
        _K = -1. * free_track_parameters::qop() * getter::norm(*_mag_field);

        // Get longitudinal momentum parallel to B field
        scalar pz = vector::dot(mom(), _h0);

        // Get transverse momentum perpendicular to B field
        vector3 pT = free_track_parameters::mom() - pz * _h0;

        // R [mm] =  pT [GeV] / B [T] in natrual unit
        _R = getter::norm(pT) / getter::norm(*_mag_field);

        // Handle the case of pT ~ 0
        if (getter::norm(pT) < 1e-6) {
            _vz_over_vt = std::numeric_limits<scalar>::infinity();
        } else {
            // Get vz over vt in new coordinate
            _vz_over_vt = pz / getter::norm(pT);
        }
    }

    /// @returns the radius of helix
    DETRAY_HOST_DEVICE
    scalar radius() const { return _R; }

    /// @returns the position after propagating the path length of s
    DETRAY_HOST_DEVICE
    point3 operator()(const scalar s) const { return this->pos(s); }

    /// @returns the position after propagating the path length of s
    DETRAY_HOST_DEVICE
    point3 pos(const scalar s) const {

        // Handle the case of pT ~ 0
        if (_vz_over_vt == std::numeric_limits<scalar>::infinity()) {
            return free_track_parameters::pos() + s * _h0;
        }

        point3 ret = free_track_parameters::pos();
        ret = ret + _delta / _K * (_K * s - std::sin(_K * s)) * _h0;
        ret = ret + std::sin(_K * s) / _K * _t0;
        ret = ret + _alpha / _K * (1 - std::cos(_K * s)) * _n0;

        return ret;
    }

    /// @returns the tangential vector after propagating the path length of s
    DETRAY_HOST_DEVICE
    vector3 dir(const scalar s) const {

        // Handle the case of pT ~ 0
        if (_vz_over_vt == std::numeric_limits<scalar>::infinity()) {
            return free_track_parameters::dir();
        }

        vector3 ret{0, 0, 0};

        ret = ret + _delta * (1 - std::cos(_K * s)) * _h0;
        ret = ret + std::cos(_K * s) * _t0;
        ret = ret + _alpha * std::sin(_K * s) * _n0;

        return ret;
    }

    /// @returns the overstep tolerance
    DETRAY_HOST_DEVICE
    scalar overstep_tolerance() const { return _overstep_tolerance; }

    /// @returns the transport jacobian after propagating the path length of s
    DETRAY_HOST_DEVICE
    free_matrix jacobian(const scalar s) const {

        free_matrix ret =
            matrix_operator().template zero<e_free_size, e_free_size>();

        const matrix_type<3, 3> I33 =
            matrix_operator().template identity<3, 3>();
        const matrix_type<3, 3> Z33 = matrix_operator().template zero<3, 3>();

        // Get drdr
        auto drdr = I33;
        matrix_operator().set_block(ret, drdr, 0, 0);

        // Get dtdr
        auto dtdr = Z33;
        matrix_operator().set_block(ret, dtdr, 4, 0);

        // Get drdt
        auto drdt = Z33;

        drdt = drdt + std::sin(_K * s) / _K * I33;

        const auto H0 = column_wise_op<scalar>().multiply(I33, _h0);
        drdt = drdt + (_K * s - std::sin(_K * s)) / _K *
                          column_wise_op<scalar>().multiply(
                              matrix_operator().transpose(H0), _h0);

        drdt = drdt + (std::cos(_K * s) - 1) / _K *
                          column_wise_op<scalar>().cross(I33, _h0);

        matrix_operator().set_block(ret, drdt, 0, 4);

        // Get dtdt
        auto dtdt = Z33;
        dtdt = dtdt + std::cos(_K * s) * I33;
        dtdt = dtdt + (1 - std::cos(_K * s)) *
                          column_wise_op<scalar>().multiply(
                              matrix_operator().transpose(H0), _h0);
        dtdt =
            dtdt - std::sin(_K * s) * column_wise_op<scalar>().cross(I33, _h0);

        matrix_operator().set_block(ret, dtdt, 4, 4);

        // Get drdl
        vector3 drdl =
            1 / free_track_parameters::qop() *
            (s * this->dir(s) + free_track_parameters::pos() - this->pos(s));

        matrix_operator().set_block(ret, drdl, 0, 7);

        // Get dtdl
        vector3 dtdl = _alpha * _K * s / free_track_parameters::qop() * _n0;

        matrix_operator().set_block(ret, dtdl, 4, 7);

        // 3x3 and 7x7 element is 1 (Maybe?)
        matrix_operator().element(ret, 3, 3) = 1;
        matrix_operator().element(ret, 7, 7) = 1;

        return ret;
    }

    private:
    /// B field
    vector3 const *_mag_field;

    /// Normalized b field
    vector3 _h0;

    /// Normalized tangent vector
    vector3 _t0;

    /// Normalized _h0 X _t0
    vector3 _n0;

    /// Magnitude of _h0 X _t0
    scalar _alpha;

    /// Dot product of _h0 X _t0
    scalar _delta;

    /// Path length scaler
    scalar _K;

    /// Radius [mm] of helix
    scalar _R;

    /// Velocity in new z axis divided by transverse velocity
    scalar _vz_over_vt;

    /// Overstep tolerance on a geomtry surface
    scalar _overstep_tolerance = -1e-4;
};

}  // namespace detail

}  // namespace detray
