/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// detray tools
#include "detray/tools/track.hpp"

namespace detray {

class helix_gun {
    public:
    using point3 = __plugin::point3<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using vector2 = __plugin::vector2<scalar>;

    helix_gun(free_track_parameters vertex, vector3 mag_field)
        : _vertex(vertex), _mag_field(mag_field) {

        // get new z direction along the b field
        _ez = vector::normalize(_mag_field);

        vector3 p = _vertex.mom();

        // get longitudinal momentum in new coordinate
        scalar pz = vector::dot(p, _ez);

        // get transverse momentum in new coordinate
        vector3 pT = p - pz * _ez;

        // R [mm] =  pT [GeV] / B [T] in natrual unit
        _R = std::abs(getter::norm(pT) / getter::norm(mag_field));

        // Handle the case of pT ~ 0
        if (getter::norm(pT) < 1e-20) {
            _vz_over_vt = std::numeric_limits<scalar>::infinity();
        } else {
            // Get vz over vt in new coordinate
            _vz_over_vt = pz / getter::norm(pT);
        }

        _scaler = _R * std::sqrt(1 + std::pow(_vz_over_vt, 2));

        // get new y direction
        scalar q = _vertex.qop() > 0 ? 1 : -1.;
        vector3 y = -1 * q * vector::normalize(pT);
        _ey = vector::normalize(y);

        // get new x direction
        _ex = vector::cross(_ey, _ez);
    }

    scalar radius() const { return _R; }

    point3 operator()(scalar path_length) const {

        // Handle the case of pT ~ 0
        if (_vz_over_vt == std::numeric_limits<scalar>::infinity()) {
            return _vertex.pos() + path_length * _ez;
        }

        // Requested path length is coverted to parameterized value (t) by
        // dividing it with scaling factor
        scalar t = path_length / _scaler;

        vector3 X = (_R * std::cos(t) - _R) * _ex;
        vector3 Y = (_R * std::sin(t)) * _ey;
        vector3 Z = (_R * _vz_over_vt * t) * _ez;

        return _vertex.pos() + X + Y + Z;
    }

    private:
    // origin of particle
    const free_track_parameters _vertex;

    // B field
    const vector3 _mag_field;

    // Radius [mm] of helix
    scalar _R;

    // new coordinated whose z axis is parallel to B field
    vector3 _ex;
    vector3 _ey;
    vector3 _ez;

    // velocity in new z axis divided by transverse velocity
    scalar _vz_over_vt;
    // The scaling factor that convert arc length into parametrized length
    scalar _scaler;
};

}  // namespace detray