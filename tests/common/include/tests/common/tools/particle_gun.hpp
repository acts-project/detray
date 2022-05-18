/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// system include
#include <cmath>
#include <utility>

// detray include(s)
#include "detray/definitions/units.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/propagator/track.hpp"
#include "detray/utils/algebra_helpers.hpp"
#include "detray/utils/enumerate.hpp"
#include "tests/common/tools/track_generators.hpp"

namespace detray {

/// @brief describes a straight-line trajectory
class ray {
    public:
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /// Parametrized constructor that complies with track interface
    ///
    /// @param pos the track position
    /// @param dir the track momentum direction
    ray(point3 pos, scalar /*time*/, vector3 dir, scalar /*q*/)
        : _pos(pos), _dir(dir) {}

    /// @returns position on the ray
    DETRAY_HOST_DEVICE point3 pos() const { return _pos; }

    /// @param position new position on the ray
    DETRAY_HOST_DEVICE void set_pos(point3 position) { _pos = position; }

    /// @returns direction of the ray
    DETRAY_HOST_DEVICE vector3 dir() const { return _dir; }

    /// @returns overstep tolerance to comply with track interface
    DETRAY_HOST_DEVICE
    scalar overstep_tolerance() const { return _overstep_tolerance; }

    private:
    /// origin of ray
    point3 _pos{0., 0., 0.};
    /// direction of ray
    vector3 _dir{0., 0., 0.};

    /// Overstep tolerance on a geomtry surface
    scalar _overstep_tolerance = -1e-4;
};

/// @brief describes a helical trajectory in a given B-field.
///
/// Helix class for the analytical solution of track propagation in
/// homogeneous B field. This Follows the notation of Eq (4.7) in
/// DOI:10.1007/978-3-030-65771-0
class helix {
    public:
    using point3 = __plugin::point3<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using vector2 = __plugin::vector2<scalar>;
    using matrix_operator = standard_matrix_operator<scalar>;
    using size_type = __plugin::size_type;
    template <size_type ROWS, size_type COLS>
    using matrix_type = __plugin::matrix_type<scalar, ROWS, COLS>;

    helix() = delete;

    /// Parametrized constructor
    ///
    /// @param vertex the helix origin
    /// @param mag_fied the magnetic field
    helix(const free_track_parameters vertex, vector3 const *const mag_field)
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
        if (getter::norm(pT) < 1e-6) {
            _vz_over_vt = std::numeric_limits<scalar>::infinity();
        } else {
            // Get vz over vt in new coordinate
            _vz_over_vt = pz / getter::norm(pT);
        }
    }

    /// @returns radius of helix
    scalar radius() const { return _R; }

    /// @returns position after propagating the path length of s
    point3 operator()(scalar s) const { return this->pos(s); }

    /// @returns position after propagating the path length of s
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

    /// @returns tangential vector after propagating the path length of s
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

    /// @returns overstep tolerance
    DETRAY_HOST_DEVICE
    scalar overstep_tolerance() const { return _overstep_tolerance; }

    /// @returns transport jacobian after propagating the path length of s
    free_matrix jacobian(scalar s) const {

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
        vector3 drdl = 1 / _vertex.qop() *
                       (s * this->dir(s) + _vertex.pos() - this->pos(s));

        matrix_operator().set_block(ret, drdl, 0, 7);

        // Get dtdl
        vector3 dtdl = _alpha * _K * s / _vertex.qop() * _n0;

        matrix_operator().set_block(ret, dtdl, 4, 7);

        // 3x3 and 7x7 element is 1 (Maybe?)
        matrix_operator().element(ret, 3, 3) = 1;
        matrix_operator().element(ret, 7, 7) = 1;

        return ret;
    }

    private:
    /// origin of particle
    free_track_parameters _vertex;

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

/// @brief struct that holds functionality to shoot a ray through a detector.
///
/// Records intersections with every detector surface along the ray.
struct particle_gun {

    using intersection_t = line_plane_intersection;

    /// Intersection implementation for straight line intersection
    template <typename surface_t, typename transform_container,
              typename mask_container>
    DETRAY_HOST_DEVICE inline static auto intersect_traj(
        const ray &r, surface_t &surface,
        const transform_container &contextual_transforms,
        const mask_container &masks) -> intersection_t {
        // Use 'intersect' function from 'intersection_kernel'
        return intersect(r, surface, contextual_transforms, masks);
    }

    /// Intersect all surfaces in a detector with a given ray.
    ///
    /// @param detector the detector
    /// @param traj the trajectory to be shot through the detector
    ///
    /// @return a sorted vector of volume indices with the corresponding
    ///         intersections of the surfaces that were encountered
    template <typename detector_t, typename trajectory_t>
    DETRAY_HOST_DEVICE inline static auto shoot_particle(
        const detector_t &detector, const trajectory_t &traj) {

        std::vector<std::pair<dindex, intersection_t>> intersection_record;

        // Loop over volumes
        for (const auto &volume : detector.volumes()) {
            for (const auto [sf_idx, sf] :
                 enumerate(detector.surfaces(), volume)) {
                // Retrieve candidate from the surface
                auto sfi = intersect_traj(traj, sf, detector.transform_store(),
                                          detector.mask_store());

                // Candidate is invalid if it oversteps too far (this is neg!)
                if (sfi.path < traj.overstep_tolerance()) {
                    continue;
                }
                // Accept if inside
                if (sfi.status == intersection::status::e_inside) {
                    // surface the candidate belongs to
                    sfi.index = volume.index();
                    intersection_record.emplace_back(sf_idx, sfi);
                }
            }
        }

        // Sort intersections by distance to origin of the trajectory
        auto sort_path = [&](std::pair<dindex, intersection_t> a,
                             std::pair<dindex, intersection_t> b) -> bool {
            return (a.second < b.second);
        };
        std::sort(intersection_record.begin(), intersection_record.end(),
                  sort_path);

        return intersection_record;
    }
};

}  // namespace detray
