/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// detray include(s)
#include "detray/propagator/track.hpp"
#include "detray/utils/algebra_helpers.hpp"

/// Helix gun class for the analytical solution of track propagation in
/// homogeneous B field. This Follows the notation of Eq (4.7) in
/// DOI:10.1007/978-3-030-65771-0
namespace detray {

class helix_gun {
    public:
    using point3 = __plugin::point3<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using vector2 = __plugin::vector2<scalar>;
    using matrix_operator = standard_matrix_operator<scalar>;
    using size_type = __plugin::size_type;
    template <size_type ROWS, size_type COLS>
    using matrix_type = __plugin::matrix_type<scalar, ROWS, COLS>;

    helix_gun() = delete;

    helix_gun(const free_track_parameters vertex,
              vector3 const *const mag_field)
        : _vertex(vertex), _mag_field(mag_field) {

        // Normalized B field
        _h0 = vector::normalize(*_mag_field);

        // Normalized tangent vector
        _t0 = vector::normalize(vertex.mom());

        // Normalized _h0 X _t0
        _n0 = vector::normalize(vector::cross(_h0, _t0));

        // Magnitude of _h0 X _t0
        _alpha = getter::norm(vector::cross(_h0, _t0));

        // Dot product of _h0 X _t0
        _delta = vector::dot(_h0, _t0);

        // Path length scaler
        _K = -1. * vertex.qop() * getter::norm(*_mag_field);

        // Get longitudinal momentum parallel to B field
        scalar pz = vector::dot(vertex.mom(), _h0);

        // Get transverse momentum perpendicular to B field
        vector3 pT = vertex.mom() - pz * _h0;

        // R [mm] =  pT [GeV] / B [T] in natrual unit
        _R = getter::norm(pT) / getter::norm(*_mag_field);

        // Handle the case of pT ~ 0
        if (getter::norm(pT) < 1e-20) {
            _vz_over_vt = std::numeric_limits<scalar>::infinity();
        } else {
            // Get vz over vt in new coordinate
            _vz_over_vt = pz / getter::norm(pT);
        }
    }

    scalar radius() const { return _R; }

    point3 operator()(scalar s) const { return this->pos(s); }

    // position after propagating the path length of s
    point3 pos(scalar s) const {

        // Handle the case of pT ~ 0
        if (_vz_over_vt == std::numeric_limits<scalar>::infinity()) {
            return _vertex.pos() + s * _h0;
        }

        point3 ret = _vertex.pos();
        ret = ret + _delta / _K * (_K * s - std::sin(_K * s)) * _h0;
        ret = ret + std::sin(_K * s) / _K * _t0;
        ret = ret + _alpha / _K * (1 - std::cos(_K * s)) * _n0;

        return ret;
    }

    // tangentional vector after propagating the path length of s
    vector3 dir(scalar s) const {

        // Handle the case of pT ~ 0
        if (_vz_over_vt == std::numeric_limits<scalar>::infinity()) {
            return _vertex.dir();
        }

        vector3 ret{0, 0, 0};

        ret = ret + _delta * (1 - std::cos(_K * s)) * _h0;
        ret = ret + std::cos(_K * s) * _t0;
        ret = ret + _alpha * std::sin(_K * s) * _n0;

        return ret;
    }

    // transport jacobian after propagating the path length of s
    free_matrix jacobian(scalar s) const {

        free_matrix ret =
            matrix_operator().template zero<e_free_size, e_free_size>();

        const matrix_type<3, 3> I33 =
            matrix_operator().template identity<3, 3>();
        const matrix_type<3, 3> Z33 = matrix_operator().template zero<3, 3>();
        /// Get drdr

        auto drdr = I33;
        matrix_operator().set_block(ret, drdr, 0, 0);

        /// Get dtdr

        auto dtdr = Z33;
        matrix_operator().set_block(ret, dtdr, 4, 0);

        /// Get drdt
        auto drdt = Z33;

        drdt = drdt + std::sin(_K * s) / _K * I33;

        const auto H0 = vector_helpers<scalar>().dot(I33, _h0);
        drdt = drdt + (_K * s - std::sin(_K * s)) / _K *
                          vector_helpers<scalar>().dot(
                              matrix_operator().transpose(H0), _h0);

        drdt = drdt + (std::cos(_K * s) - 1) / _K *
                          vector_helpers<scalar>().cross(I33, _h0);

        matrix_operator().set_block(ret, drdt, 0, 4);

        /// Get dtdt
        auto dtdt = Z33;
        dtdt = dtdt + std::cos(_K * s) * I33;
        dtdt = dtdt + (1 - std::cos(_K * s)) *
                          vector_helpers<scalar>().dot(
                              matrix_operator().transpose(H0), _h0);
        dtdt =
            dtdt - std::sin(_K * s) * vector_helpers<scalar>().cross(I33, _h0);

        matrix_operator().set_block(ret, dtdt, 4, 4);

        /// Get drdl
        vector3 drdl = 1 / _vertex.qop() *
                       (s * this->dir(s) + _vertex.pos() - this->pos(s));

        matrix_operator().set_block(ret, drdl, 0, 7);

        /// Get dtdl
        vector3 dtdl = _alpha * _K * s / _vertex.qop() * _n0;

        matrix_operator().set_block(ret, dtdl, 4, 7);

        /// 3x3 and 7x7 element is 1 (Maybe?)
        matrix_operator().element(ret, 3, 3) = 1;
        matrix_operator().element(ret, 7, 7) = 1;

        return ret;
    }

    private:
    // origin of particle
    free_track_parameters _vertex;

    // B field
    vector3 const *_mag_field;

    // Normalized b field
    vector3 _h0;

    // Normalized tangent vector
    vector3 _t0;

    // Normalized _h0 X _t0
    vector3 _n0;

    // Magnitude of _h0 X _t0
    scalar _alpha;

    // Dot product of _h0 X _t0
    scalar _delta;

    // Path length scaler
    scalar _K;

    // Radius [mm] of helix
    scalar _R;

    // Velocity in new z axis divided by transverse velocity
    scalar _vz_over_vt;
};

}  // namespace detray
