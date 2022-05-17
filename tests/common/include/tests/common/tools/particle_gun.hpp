/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cmath>
#include <utility>

#include "detray/definitions/units.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/utils/enumerate.hpp"

namespace detray {

/// Track parametrization that reflects a simple straight-line track
struct ray {
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /// Generate default constructor
    ray() = default;

    /// Parametrized constructor that complies with track interface
    ///
    /// @param pos the track position
    /// @param dir the track momentum direction
    ray(point3 pos, scalar /*time*/, vector3 dir, scalar /*q*/)
        : _pos(pos), _dir(dir) {}

    /// Current position in the geometry during propagation
    point3 _pos = {0., 0., 0.};
    /// Diretion
    vector3 _dir = {0., 0., 0.};
    /// Overstep tolerance on a geomtry surface (not needed here)
    scalar _overstep_tolerance = -1e-4;

    /// @returns position on the ray
    DETRAY_HOST_DEVICE point3 pos() const { return _pos; }

    /// @returns direction of the ray
    DETRAY_HOST_DEVICE vector3 dir() const { return _dir; }

    /// @returns overstep tolerance to comply with track interface
    DETRAY_HOST_DEVICE
    scalar overstep_tolerance() const { return _overstep_tolerance; }
};

/// @brief Genrates track states with momentum directions in a uniform angle
/// space.
///
/// It generates the track instances on the fly according to given parameters
/// and with the momentum direction determined by phi and theta angles, which
/// are generated as the iteration proceeds. The stepsizes in the angle space
/// spans theta ]0,pi[ x phi [-pi, pi], while the step sizes (and with them
/// the number of generated tracks) are configurable.
///
/// @tparam track_t the type of track parametrization that should be used.
template <typename track_t>
struct uniform_track_generator {
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /// Start and end of angle space
    std::size_t m_theta_steps{50};
    std::size_t m_phi_steps{50};
    scalar m_phi{0}, m_theta{0};

    /// Track params
    point3 m_origin{};
    scalar m_mom_mag = 10. * unit_constants::GeV;

    /// Iteration indices
    std::size_t i_phi{0}, i_theta{0};

    /// Default constructor
    uniform_track_generator() = default;

    /// Paramtetrized constructor for fine-grained configurations
    ///
    /// @param theta_steps the number of steps in the theta space
    /// @param phi_steps the number of steps in the phi space
    /// @param trk_origin the starting point of the track
    /// @param trk_mom magnitude of the track momentum
    uniform_track_generator(std::size_t theta_steps, std::size_t phi_steps,
                            point3 trk_origin = {},
                            scalar trk_mom = 10. * unit_constants::GeV)
        : m_theta_steps(theta_steps),
          m_phi_steps(phi_steps),
          m_origin(trk_origin),
          m_mom_mag(trk_mom),
          i_phi(0),
          i_theta(0) {}

    /// @returns generator in starting state
    DETRAY_HOST_DEVICE
    auto begin() {
        i_phi = 0;
        i_theta = 0;
        return *this;
    }

    /// @returns generator in end state
    DETRAY_HOST_DEVICE
    auto end() {
        i_phi = m_phi_steps;
        i_theta = m_theta_steps;
        return *this;
    }

    /// @returns whether we reached end of angle space
    DETRAY_HOST_DEVICE
    bool operator!=(const uniform_track_generator &rhs) const {
        return not(rhs.m_theta_steps == m_theta_steps and
                   rhs.m_phi_steps == m_phi_steps and rhs.i_phi == i_phi and
                   rhs.i_theta == i_theta);
    }

    /// Iterate through angle space according to given step sizes.
    DETRAY_HOST_DEVICE
    void operator++() {
        // Check theta range according to step size
        if (i_theta < m_theta_steps) {
            // Calculate theta ]0,pi[
            m_theta = 0.01 + i_theta * (M_PI - 0.01) / m_theta_steps;
            // Check phi sub-range
            if (i_phi < m_phi_steps) {
                // Calculate phi [-pi, pi]
                m_phi = -M_PI + i_phi * (2 * M_PI) / m_phi_steps;
                ++i_phi;
                // Don't update theta while phi sub-range is not done
                return;
            }
            // Reset phi range
            i_phi = 0;
            ++i_theta;
        } else {
            // This got reset to zero in last pass
            i_phi = m_phi_steps;
        }
    }

    /// Genrate the track instance
    DETRAY_HOST_DEVICE
    track_t operator*() const {
        // Momentum direction from angles
        vector3 mom{std::cos(m_phi) * std::sin(m_theta),
                    std::sin(m_phi) * std::sin(m_theta), std::cos(m_theta)};
        // Magnitude of momentum
        vector::normalize(mom);
        mom = m_mom_mag * mom;

        return track_t{m_origin, 0, mom, -1};
    }
};

/// Intersect all portals in a detector with a given ray.
///
/// @param detector the detector
/// @param ray the ray to be shot through the detector
///
/// @return a sorted vector of volume indices with the corresponding
///         intersections of the portals that were encountered
template <typename detector_t>
inline auto shoot_ray(const detector_t &detector, const ray &r) {

    using intersection_t = line_plane_intersection;

    std::vector<std::pair<dindex, intersection_t>> intersection_record;

    // Loop over volumes
    for (const auto &volume : detector.volumes()) {
        for (const auto [sf_idx, sf] : enumerate(detector.surfaces(), volume)) {
            // Retrieve candidate from the object
            auto sfi = intersect(r, sf, detector.transform_store(),
                                 detector.mask_store());

            // Candidate is invalid if it oversteps too far (this is neg!)
            if (sfi.path < r.overstep_tolerance()) {
                continue;
            }
            // Accept if inside
            if (sfi.status == intersection::status::e_inside) {
                // object the candidate belongs to
                sfi.index = volume.index();
                intersection_record.emplace_back(sf_idx, sfi);
            }
        }
    }

    // Sort intersections by distance to origin of the ray
    auto sort_path = [&](std::pair<dindex, intersection_t> a,
                         std::pair<dindex, intersection_t> b) -> bool {
        return (a.second < b.second);
    };
    std::sort(intersection_record.begin(), intersection_record.end(),
              sort_path);

    return intersection_record;
}

/// Intersect all portals in a detector with a given ray.
///
/// @param detector the detector
/// @param origin the origin of the ray in global coordinates
/// @param direction the direction of the ray in global coordinater
///
/// @return a sorted vector of volume indices with the corresponding
///         intersections of the portals that were encountered
template <typename detector_t, typename point_t>
inline auto shoot_ray(const detector_t &detector, const point_t &origin,
                      const point_t &direction) {
    ray r{origin, direction};
    return shoot_ray(detector, r);
}

}  // namespace detray
