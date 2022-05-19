/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cmath>
#include <iostream>
#include <utility>

#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"

namespace detray {

/// @brief Generates track states with momentum directions in a uniform angle
/// space.
///
/// It generates the track instances on the fly according to given parameters
/// and with the momentum direction determined by phi and theta angles, which
/// are advanced as the iteration proceeds. The angle space spans
/// theta ]0, pi[ x phi [-pi, pi], while the step sizes (and with them
/// the number of generated tracks) are configurable.
///
/// @tparam track_t the type of track parametrization that should be used.
template <typename track_t>
class uniform_track_generator {
    public:
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /// Default constructor
    DETRAY_HOST_DEVICE
    uniform_track_generator() = default;

    /// Paramtetrized constructor for fine-grained configurations
    ///
    /// @param theta_steps the number of steps in the theta space
    /// @param phi_steps the number of steps in the phi space
    /// @param trk_origin the starting point of the track
    /// @param trk_mom magnitude of the track momentum
    DETRAY_HOST_DEVICE
    uniform_track_generator(std::size_t n_theta, std::size_t n_phi,
                            point3 trk_origin = {},
                            scalar trk_mom = 1. * unit_constants::GeV,
                            scalar time = 0., scalar charge = -1.)
        : m_theta_steps(n_theta),
          m_phi_steps(n_phi),
          m_origin(trk_origin),
          m_mom_mag(trk_mom),
          m_time{time},
          m_charge{charge},
          i_phi(0),
          i_theta(0) {}

    /// @returns generator in starting state: Default values resolve the first
    /// phi angle iteration.
    DETRAY_HOST_DEVICE
    auto begin() -> uniform_track_generator {
        i_phi = 1;
        i_theta = 0;
        return *this;
    }

    /// @returns generator in end state
    DETRAY_HOST_DEVICE
    auto end() -> uniform_track_generator {
        i_phi = 1;
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
    auto operator++() -> uniform_track_generator {
        scalar pi{M_PI};
        // Check theta range according to step size
        if (i_theta < m_theta_steps) {
            // Check phi sub-range
            if (i_phi < m_phi_steps) {
                // Calculate new phi [-pi, pi]
                m_phi = -pi + i_phi * (scalar{2} * pi) / m_phi_steps;
                ++i_phi;
                return *this;
            }
            // Reset phi range
            i_phi = 1;
            m_phi = -M_PI;
            // Calculate new theta ]0,pi[
            ++i_theta;
            m_theta =
                scalar{0.01} + i_theta * (pi - scalar{0.01}) / m_theta_steps;
        }
        return *this;
    }

    /// @return track instance from generated momentum direction
    DETRAY_HOST_DEVICE
    track_t operator*() const {
        // Momentum direction from angles
        vector3 mom{std::cos(m_phi) * std::sin(m_theta),
                    std::sin(m_phi) * std::sin(m_theta), std::cos(m_theta)};
        // Magnitude of momentum
        vector::normalize(mom);
        mom = m_mom_mag * mom;

        return track_t{m_origin, m_time, mom, m_charge};
    }

    protected:
    /// Start and end of angle space
    std::size_t m_theta_steps{50};
    std::size_t m_phi_steps{50};
    /// Phi and theta angles of momentum direction
    scalar m_phi{-M_PI}, m_theta{0.01};

    /// Track origin
    point3 m_origin{0., 0., 0.};
    /// Magnitude of momentum: Default is one to keep directions normalized if
    /// no momentum information is needed (e.g. for a ray)
    scalar m_mom_mag = 1. * unit_constants::GeV;
    /// Time parameter and charge of the track
    scalar m_time{0}, m_charge{0};

    /// Iteration indices
    std::size_t i_phi{0};
    std::size_t i_theta{0};
};

/// @brief Generates track states with random momentum directions.
///
///
///
/// @tparam track_t the type of track parametrization that should be used.
/*template <typename track_t, typename random_nr_generator_t>
struct random_track_generator {
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /// Start and end of angle space
    scalar m_phi{0}, m_theta{0};

    /// Track params
    point3 m_origin{0., 0., 0.};
    /// Magnitude of momentum: Default is one to keep directions normalized if
    /// no momentum information is needed (e.g. for a ray)
    scalar m_mom_mag = 1. * unit_constants::GeV;

    /// Iteration indices
    std::size_t m_n_tracks{50*50};
    std::size_t i_track{0};

    /// Default constructor
    DETRAY_HOST_DEVICE
    random_track_generator() = default;

    /// Paramtetrized constructor for fine-grained configurations
    ///
    /// @param theta_steps the number of steps in the theta space
    /// @param phi_steps the number of steps in the phi space
    /// @param trk_origin the starting point of the track
    /// @param trk_mom magnitude of the track momentum
    DETRAY_HOST_DEVICE
    random_track_generator(std::size_t n_tracks, point3 trk_origin = {},
                           scalar trk_mom = 1. * unit_constants::GeV)
        : m_n_tracks(n_tracks),
          m_origin(trk_origin),
          m_mom_mag(trk_mom),
          i_track(0)
          {}

    /// @returns generator in starting state
    DETRAY_HOST_DEVICE
    auto begin() {
        i_track = 0;
        return *this;
    }

    /// @returns generator in end state
    DETRAY_HOST_DEVICE
    auto end() {
        i_track = m_n_tracks;
        return *this;
    }

    /// @returns whether we reached end of angle space
    DETRAY_HOST_DEVICE
    bool operator!=(const uniform_track_generator &rhs) const {
        return not(rhs.m_n_tracks == m_n_tracks and rhs.i_track == i_track);
    }

    /// Iterate through angle space according to given step sizes.
    DETRAY_HOST_DEVICE
    void operator++() {
        ...
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
};*/

}  // namespace detray
