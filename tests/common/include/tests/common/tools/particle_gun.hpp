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
class helix : free_track_parameters {
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
    /// @param pos the the origin of the helix
    /// @param time the time parameter
    /// @param dir the initial direction of momentum for the helix
    /// @param q the charge of the particle
    /// @param mag_field the magnetic field vector
    helix(point3 pos, scalar time, vector3 dir, scalar q,
          vector3 const *const mag_field)
        : free_track_parameters(pos, time, dir, q), _mag_field(mag_field) {}

    /// Parametrized constructor
    ///
    /// @param vertex the underlying track parametrization
    /// @param mag_fied the magnetic field vector
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

    /// @returns radius of helix
    scalar radius() const { return _R; }

    /// @returns position after propagating the path length of s
    point3 operator()(scalar s) const { return this->pos(s); }

    /// @returns position after propagating the path length of s
    point3 pos(scalar s) const {

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

    /// @returns tangential vector after propagating the path length of s
    vector3 dir(scalar s) const {

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

/// @brief struct that holds functionality to shoot a ray through a detector.
///
/// Records intersections with every detector surface along the ray.
struct particle_gun {

    using intersection_t = line_plane_intersection;

    /// Intersect all surfaces in a detector with a given ray.
    ///
    /// @param detector the detector.
    /// @param traj the trajectory to be shot through the detector.
    ///
    /// @return a sorted vector of volume indices with the corresponding
    ///         intersections of the surfaces that were encountered.
    template <typename detector_t, typename trajectory_t>
    DETRAY_HOST_DEVICE inline static auto shoot_particle(
        const detector_t &detector, const trajectory_t &traj) {

        std::vector<std::pair<dindex, intersection_t>> intersection_record;

        // Loop over volumes
        for (const auto &volume : detector.volumes()) {
            for (const auto [sf_idx, sf] :
                 enumerate(detector.surfaces(), volume)) {
                // Retrieve candidate from the surface
                auto sfi = intersect(traj, sf, detector.transform_store(),
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

    /// Wrap the @c intersection_kernel intersect function for rays
    template <typename surface_t, typename transform_container,
              typename mask_container>
    DETRAY_HOST_DEVICE static inline auto intersect(
        const ray &r, surface_t &surface,
        const transform_container &contextual_transforms,
        const mask_container &masks) {
        return detray::intersect(r, surface, contextual_transforms, masks);
    }

    /// Start helix intersection unrolling. See documentation of the detray
    /// intersection kernel
    template <typename surface_t, typename transform_container,
              typename mask_container>
    DETRAY_HOST_DEVICE static inline auto intersect(
        const helix &h, surface_t &surface,
        const transform_container &contextual_transforms,
        const mask_container &masks) {
        // Gather all information to perform intersections
        const auto &ctf = contextual_transforms[surface.transform()];
        const auto volume_index = surface.volume();
        const auto mask_id = surface.mask_type();
        const auto &mask_range = surface.mask_range();

        // No intersectors available for non-planar surfaces
        if (surface.is_portal()) {
            return intersection_t{};
        }
        // Unroll the intersection depending on the mask container size
        using mask_defs = typename surface_t::mask_defs;

        return unroll_helix_intersect(
            h, ctf, masks, mask_range, mask_id, volume_index,
            std::make_integer_sequence<unsigned int, mask_defs::n_types>{});
    }

    /// Helix version of the intersection unrolling. Calls the
    /// @c helix_intersector of this class instead of the mask's native
    /// intersector.
    template <typename transform_t, typename mask_container_t,
              typename mask_range_t, unsigned int first_mask_id,
              unsigned int... remaining_mask_ids>
    DETRAY_HOST_DEVICE static inline auto unroll_helix_intersect(
        const helix &h, const transform_t &ctf, const mask_container_t &masks,
        const mask_range_t &mask_range,
        const typename mask_container_t::id_type mask_id, dindex volume_index,
        std::integer_sequence<unsigned int, first_mask_id,
                              remaining_mask_ids...>
        /*available_ids*/) {

        // Pick the first one for interseciton
        if (mask_id == first_mask_id) {

            auto &mask_group =
                masks.template group<mask_container_t::to_id(first_mask_id)>();

            // Check all masks of this surface for intersection with the helix
            for (const auto &mask : range(mask_group, mask_range)) {
                auto sfi = std::move(helix_intersect(ctf, h, mask));

                if (sfi.status == intersection::status::e_inside) {
                    sfi.index = volume_index;
                    return sfi;
                }
            }
        }

        // The reduced integer sequence
        std::integer_sequence<unsigned int, remaining_mask_ids...> remaining;

        // Unroll as long as you have at least 1 entries
        if constexpr (remaining.size() >= 1) {
            return (unroll_helix_intersect(h, ctf, masks, mask_range, mask_id,
                                           volume_index, remaining));
        }

        // No intersection was found
        return intersection_t{};
    }

    /// @brief Intersection implementation for helical trajectories.
    ///
    /// The algorithm uses the Newton-Raphson method to find an intersection on
    /// the unbounded surface and then applies the mask.
    ///
    /// @return the intersection.
    template <typename transform_t, typename mask_t>
    DETRAY_HOST_DEVICE inline static auto helix_intersector(
        const transform_t &trf, const helix &h, const mask_t &mask,
        const typename mask_t::mask_tolerance tolerance =
            mask_t::within_epsilon) -> intersection_t {

        using local_frame = typename mask_t::local_type;

        // Get the surface info
        const auto &sm = trf.matrix();
        auto sn = getter::vector<3>(sm, 0, 2);
        auto st = getter::vector<3>(sm, 0, 3);

        // starting point on the helix for the Newton iteration
        scalar epsilon = 0.001;
        scalar s{getter::norm(st)};

        // Guard against inifinite loops
        std::size_t n_tries{0};
        std::size_t max_n_tries = 1000;

        // Run the iteration on s
        while (s > epsilon and n_tries < max_n_tries) {
            s -= vector::dot(sn, h.pos(s) - st) / vector::dot(sn, h.dir(s));
            ++n_tries;
        }
        // No intersection found within max number of trials
        if (n_tries == max_n_tries) {
            return intersection_t{};
        }

        // Build intersection struct from helix parameter s
        intersection_t is;
        is.path = getter::norm(h.pos(s));
        is.p3 = h.pos(s);
        constexpr local_frame local_converter{};
        is.p2 = local_converter(trf, is.p3);
        is.status = mask.template is_inside<local_frame>(is.p2, tolerance);
        is.direction = is.path > h.overstep_tolerance()
                           ? intersection::direction::e_along
                           : intersection::direction::e_opposite;
        is.link = mask.volume_link();
        return is;
    }
};

}  // namespace detray
