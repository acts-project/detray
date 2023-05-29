/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/tracks/detail/track_helper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/matrix_helper.hpp"

// System include(s).
#include <cmath>
#include <limits>
#include <ostream>

namespace detray {

namespace detail {

/// @brief describes a straight-line trajectory
template <typename transform3_t>
class simple_ray {
    public:
    using transform3_type = transform3_t;
    using scalar_type = typename transform3_type::scalar_type;
    using vector3 = typename transform3_type::vector3;
    using point3 = typename transform3_type::point3;

    /// Parametrized constructor that complies with track interface
    ///
    /// @param pos the track position
    /// @param dir the track momentum direction
    DETRAY_HOST_DEVICE simple_ray(const point3 &pos, const vector3 &dir,
                                  const scalar_type /*t*/,
                                  const scalar_type /*qop*/)
        : _pos{pos}, _dir{dir} {}

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

    /// Print
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &os, const simple_ray &r) {
        os << "ray: ";
        os << "ori = [" << r._pos[0] << ", " << r._pos[1] << ", " << r._pos[2]
           << "], ";
        os << "dir = [" << r._dir[0] << ", " << r._dir[1] << ", " << r._dir[2]
           << "], ";
        os << "tol = " << r._overstep_tolerance << std::endl;

        return os;
    }

    private:
    /// origin of ray
    point3 _pos{0.f, 0.f, 0.f};
    /// direction of ray
    vector3 _dir{0.f, 0.f, 1.f};

    /// Overstep tolerance on a geometry surface
    scalar_type _overstep_tolerance{-1e-4f};
};

/// @brief describes a straight-line trajectory
template <typename transform3_t>
class ray {
    public:
    using transform3_type = transform3_t;
    using scalar_type = typename transform3_type::scalar_type;
    using vector3 = typename transform3_type::vector3;
    using point3 = typename transform3_type::point3;

    using free_track_parameters_type = free_track_parameters<transform3_type>;
    using free_vector_type = typename free_track_parameters_type::vector_type;

    // Track helper
    using matrix_operator = typename transform3_t::matrix_actor;
    using track_helper = detail::track_helper<matrix_operator>;

    ray() = default;

    /// Parametrized constructor that complies with track interface
    ///
    /// @param pos the track position
    /// @param dir the track momentum direction
    DETRAY_HOST_DEVICE ray(const point3 &pos, const scalar_type /*time*/,
                           const vector3 &dir, const scalar_type /*qop*/)
        : _pos{pos}, _dir{dir} {}

    /// Parametrized constructor that complies with track interface
    ///
    /// @param track the track state that should be approximated
    DETRAY_HOST_DEVICE ray(const free_vector_type &free_vec)
        : ray(track_helper().pos(free_vec), track_helper().time(free_vec),
              track_helper().dir(free_vec), track_helper().qop(free_vec)) {}

    /// Parametrized constructor that complies with track interface
    ///
    /// @param track the track state that should be approximated
    DETRAY_HOST_DEVICE ray(const free_track_parameters_type &track)
        : ray(track.vector()) {
        _overstep_tolerance = track.overstep_tolerance();
    }

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

    /// Print
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &os, const ray &r) {
        os << "ray: ";
        os << "ori = [" << r._pos[0] << ", " << r._pos[1] << ", " << r._pos[2]
           << "], ";
        os << "dir = [" << r._dir[0] << ", " << r._dir[1] << ", " << r._dir[2]
           << "], ";
        os << "tol = " << r._overstep_tolerance << std::endl;

        return os;
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
class helix {
    public:
    using transform3_type = transform3_t;
    using scalar_type = typename transform3_type::scalar_type;
    using matrix_operator = typename transform3_type::matrix_actor;
    using vector3 = typename transform3_type::vector3;
    using point3 = typename transform3_type::point3;

    /// Free track parameters
    using free_track_parameters_type = free_track_parameters<transform3_t>;
    using free_vector_type = typename free_track_parameters_type::vector_type;

    /// Size type
    using size_type = typename transform3_type::size_type;
    /// 2D Matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    using free_matrix = matrix_type<e_free_size, e_free_size>;
    using mat_helper = matrix_helper<matrix_operator>;

    // Track helper
    using track_helper = detail::track_helper<matrix_operator>;

    DETRAY_HOST_DEVICE
    helix() = delete;

    /// Parametrized constructor
    ///
    /// @param pos the the origin of the helix
    /// @param time the time parameter
    /// @param dir the initial direction for the helix
    /// @param q the charge of the particle
    /// @param mag_field the magnetic field vector
    DETRAY_HOST_DEVICE
    helix(const point3 &pos, const scalar_type time, const vector3 &dir,
          const scalar_type qop, vector3 const *const mag_field)
        : _pos(pos), _time(time), _qop(qop), _mag_field(mag_field) {

        // Normalized B field
        _h0 = vector::normalize(*_mag_field);

        // Normalized tangent vector
        _t0 = dir;

        // Momentum
        const vector3 mom = 1.f / static_cast<scalar_type>(std::abs(qop)) * dir;

        // Normalized _h0 X _t0
        _n0 = vector::normalize(vector::cross(_h0, _t0));

        // Magnitude of _h0 X _t0
        _alpha = getter::norm(vector::cross(_h0, _t0));

        // Dot product of _h0 X _t0
        _delta = vector::dot(_h0, _t0);

        // Path length scaler
        _K = -_qop * getter::norm(*_mag_field);

        // Get longitudinal momentum parallel to B field
        scalar_type pz = vector::dot(mom, _h0);

        // Get transverse momentum perpendicular to B field
        vector3 pT = mom - pz * _h0;

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

    DETRAY_HOST_DEVICE
    helix(const free_vector_type &free_vec, vector3 const *const mag_field)
        : helix(track_helper().pos(free_vec), track_helper().time(free_vec),
                track_helper().dir(free_vec), track_helper().qop(free_vec),
                mag_field) {}

    DETRAY_HOST_DEVICE
    helix(const free_track_parameters_type &track,
          vector3 const *const mag_field)
        : helix(track.vector(), mag_field) {}

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
            return _pos + s * _h0;
        }

        point3 ret = _pos;
        ret = ret + _delta / _K * (_K * s - math_ns::sin(_K * s)) * _h0;
        ret = ret + math_ns::sin(_K * s) / _K * _t0;
        ret = ret + _alpha / _K * (1.f - math_ns::cos(_K * s)) * _n0;

        return ret;
    }

    DETRAY_HOST_DEVICE
    point3 pos() const { return _pos; }

    /// @returns the tangential vector after propagating the path length of s
    DETRAY_HOST_DEVICE
    vector3 dir(const scalar_type s) const {

        // Handle the case of pT ~ 0
        if (_vz_over_vt == std::numeric_limits<scalar_type>::infinity()) {
            return _t0;
        }

        vector3 ret{0.f, 0.f, 0.f};

        ret = ret + _delta * (1 - math_ns::cos(_K * s)) * _h0;
        ret = ret + math_ns::cos(_K * s) * _t0;
        ret = ret + _alpha * math_ns::sin(_K * s) * _n0;

        return ret;
    }

    DETRAY_HOST_DEVICE
    point3 dir() const { return _t0; }

    DETRAY_HOST_DEVICE
    scalar_type time() const { return _time; }

    DETRAY_HOST_DEVICE
    scalar_type qop() const { return _qop; }

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

        drdt = drdt + math_ns::sin(_K * s) / _K * I33;

        const auto H0 = mat_helper().column_wise_multiply(I33, _h0);
        drdt = drdt + (_K * s - math_ns::sin(_K * s)) / _K *
                          mat_helper().column_wise_multiply(
                              matrix_operator().transpose(H0), _h0);

        drdt = drdt + (math_ns::cos(_K * s) - 1.f) / _K *
                          mat_helper().column_wise_cross(I33, _h0);

        matrix_operator().set_block(ret, drdt, e_free_pos0, e_free_dir0);

        // Get dtdt
        auto dtdt = Z33;
        dtdt = dtdt + math_ns::cos(_K * s) * I33;
        dtdt = dtdt + (1 - math_ns::cos(_K * s)) *
                          mat_helper().column_wise_multiply(
                              matrix_operator().transpose(H0), _h0);
        dtdt = dtdt -
               math_ns::sin(_K * s) * mat_helper().column_wise_cross(I33, _h0);

        matrix_operator().set_block(ret, dtdt, e_free_dir0, e_free_dir0);

        // Get drdl
        vector3 drdl = 1.f / _qop * (s * this->dir(s) + _pos - this->pos(s));

        matrix_operator().set_block(ret, drdl, e_free_pos0, e_free_qoverp);

        // Get dtdl
        vector3 dtdl = _alpha * _K * s / _qop * _n0;

        matrix_operator().set_block(ret, dtdl, e_free_dir0, e_free_qoverp);

        // 3x3 and 7x7 element is 1 (Maybe?)
        matrix_operator().element(ret, e_free_time, e_free_time) = 1.f;
        matrix_operator().element(ret, e_free_qoverp, e_free_qoverp) = 1.f;

        return ret;
    }

    /// Print
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &os, const helix &h) {
        os << "helix: ";
        os << "ori = [" << h._pos[0] << ", " << h._pos[1] << ", " << h._pos[2]
           << "], ";
        os << "dir = [" << h._t0[0] << ", " << h._t0[1] << ", " << h._t0[2]
           << "]" << std::endl;

        return os;
    }

    /// origin
    point3 _pos;

    /// time
    scalar_type _time;

    /// qop
    scalar_type _qop;

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

}  // namespace detail

template <typename algebra_t>
using ray = detail::simple_ray<algebra_t>;

}  // namespace detray
