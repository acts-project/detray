/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/concepts/track.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/utils/ranges/ranges.hpp"

// System include(s)
#include <algorithm>
#include <cmath>

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
template <CONSTRAINT(concepts::track) track_t>
class uniform_track_generator
    : public detray::ranges::view_interface<uniform_track_generator<track_t>> {
    public:
    using point3 = typename track_t::point3;
    using vector3 = typename track_t::vector3;

    /// Configure how tracks are generated
    struct configuration {

        /// Range for theta and phi
        std::array<scalar, 2> m_theta_range{0.01f, constant<scalar>::pi};
        std::array<scalar, 2> m_phi_range{-constant<scalar>::pi,
                                          constant<scalar>::pi};
        /// Angular step size
        std::size_t m_theta_steps{50u};
        std::size_t m_phi_steps{50u};

        /// Track origin
        point3 m_origin{0.f, 0.f, 0.f};

        /// Magnitude of momentum: Default is one to keep directions normalized
        /// if no momentum information is needed (e.g. for a ray)
        scalar m_p_mag{1.f * unit<scalar>::GeV};

        /// Time parameter and charge of the track
        scalar m_time{0.f * unit<scalar>::us}, m_charge{-1.f * unit<scalar>::e};

        /// Setters
        /// @{
        configuration& theta_range(scalar low, scalar high) {
            m_theta_range = {low, high};
            return *this;
        }
        configuration& phi_range(scalar low, scalar high) {
            m_phi_range = {low, high};
            return *this;
        }
        configuration& theta_steps(std::size_t n) {
            m_theta_steps = n;
            return *this;
        }
        configuration& phi_steps(std::size_t n) {
            m_phi_steps = n;
            return *this;
        }
        configuration& origin(point3 ori) {
            m_origin = ori;
            return *this;
        }
        configuration& p_mag(scalar p) {
            m_p_mag = p;
            return *this;
        }
        configuration& time(scalar t) {
            m_time = t;
            return *this;
        }
        configuration& charge(scalar q) {
            m_charge = q;
            return *this;
        }
        /// @}

        /// Getters
        /// @{
        constexpr std::array<scalar, 2> theta_range() const {
            return m_theta_range;
        }
        constexpr std::array<scalar, 2> phi_range() const {
            return m_phi_range;
        }
        constexpr std::size_t theta_steps() const { return m_theta_steps; }
        constexpr std::size_t phi_steps() const { return m_phi_steps; }
        constexpr const point3& origin() const { return m_origin; }
        constexpr scalar p_mag() const { return m_p_mag; }
        constexpr scalar time() const { return m_time; }
        constexpr scalar charge() const { return m_charge; }
        /// @}
    };

    private:
    /// @brief Nested iterator type that generates track states.
    struct iterator {

        using difference_type = std::ptrdiff_t;
        using value_type = track_t;
        using pointer = track_t*;
        using reference = track_t&;
        using iterator_category = detray::ranges::input_iterator_tag;

        constexpr iterator() = default;

        DETRAY_HOST_DEVICE
        constexpr iterator(configuration cfg, std::size_t iph = 1u,
                           std::size_t ith = 0u)
            : m_cfg{cfg},
              m_theta_step_size{(cfg.theta_range()[1] - cfg.theta_range()[0]) /
                                static_cast<scalar>(cfg.theta_steps())},
              m_phi_step_size{(cfg.phi_range()[1] - cfg.phi_range()[0]) /
                              static_cast<scalar>(cfg.phi_steps())},
              m_phi{cfg.phi_range()[0]},
              m_theta{cfg.theta_range()[0]},
              i_phi{iph},
              i_theta{ith} {}

        /// @returns whether we reached end of angle space
        DETRAY_HOST_DEVICE
        constexpr bool operator==(const iterator& rhs) const {
            return rhs.i_phi == i_phi and rhs.i_theta == i_theta;
        }

        /// @returns whether we reached end of angle space
        DETRAY_HOST_DEVICE
        constexpr bool operator!=(const iterator& rhs) const {
            return not(*this == rhs);
        }

        /// Iterate through angle space according to given step sizes.
        ///
        /// @returns the generator at its next position.
        DETRAY_HOST_DEVICE
        constexpr auto operator++() -> iterator& {
            // Check theta range according to step size
            if (i_theta < m_cfg.theta_steps()) {
                // Check phi sub-range
                if (i_phi < m_cfg.phi_steps()) {
                    // Calculate new phi in the given range
                    m_phi = m_cfg.phi_range()[0] +
                            static_cast<scalar>(i_phi) * m_phi_step_size;
                    ++i_phi;
                    return *this;
                }
                // Reset phi range
                i_phi = 1;
                m_phi = m_cfg.phi_range()[0];

                // Calculate new theta in the given range
                ++i_theta;
                m_theta = m_cfg.theta_range()[0] +
                          static_cast<scalar>(i_theta) * m_theta_step_size;
            }
            return *this;
        }

        /// @returns a track instance from generated momentum direction
        DETRAY_HOST_DEVICE
        track_t operator*() const {
            // Momentum direction from angles
            vector3 mom{math_ns::cos(m_phi) * std::sin(m_theta),
                        std::sin(m_phi) * std::sin(m_theta),
                        math_ns::cos(m_theta)};
            // Magnitude of momentum
            vector::normalize(mom);
            mom = m_cfg.p_mag() * mom;

            return track_t{m_cfg.origin(), m_cfg.time(), mom, m_cfg.charge()};
        }

        /// Current configuration
        configuration m_cfg{};

        /// Angular step sizes
        scalar m_theta_step_size{0.f};
        scalar m_phi_step_size{0.f};

        /// Phi and theta angles of momentum direction
        scalar m_phi{-constant<scalar>::pi}, m_theta{0.01f};

        /// Iteration indices
        std::size_t i_phi{0u};
        std::size_t i_theta{0u};
    };

    configuration m_cfg{};

    public:
    using iterator_t = iterator;

    /// Default constructor
    constexpr uniform_track_generator() = default;

    /// Construct from external configuration @param cfg
    DETRAY_HOST_DEVICE
    constexpr uniform_track_generator(configuration cfg) : m_cfg{cfg} {}

    /// Paramtetrized constructor for fine-grained configurations
    ///
    /// @param n_theta the number of steps in the theta space
    /// @param n_phi the number of steps in the phi space
    /// @param trk_origin the starting point of the track
    /// @param trk_mom magnitude of the track momentum (in GeV)
    /// @param theta_range the range for theta values
    /// @param phi_range the range for phi values
    /// @param time time measurement (micro seconds)
    /// @param charge charge of particle (e)
    DETRAY_HOST_DEVICE
    uniform_track_generator(
        std::size_t n_theta, std::size_t n_phi,
        point3 trk_origin = {0.f, 0.f, 0.f},
        scalar trk_mom = 1.f * unit<scalar>::GeV,
        std::array<scalar, 2> theta_range = {0.01f, constant<scalar>::pi},
        std::array<scalar, 2> phi_range = {-constant<scalar>::pi,
                                           constant<scalar>::pi},
        scalar time = 0.f * unit<scalar>::us,
        scalar charge = -1.f * unit<scalar>::e)
        : m_cfg{theta_range, phi_range, n_theta, n_phi,
                trk_origin,  trk_mom,   time,    charge} {}

    /// Move constructor
    uniform_track_generator(uniform_track_generator&& other)
        : m_cfg(std::move(other.m_cfg)) {}

    /// Copy assignment operator
    DETRAY_HOST_DEVICE
    constexpr uniform_track_generator& operator=(
        const uniform_track_generator& other) {
        m_cfg = other.m_cfg;
        return *this;
    }

    /// Access the configuration
    constexpr configuration& config() { return m_cfg; }

    /// @returns the generator in initial state: Default values reflect the
    /// first phi angle iteration.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const noexcept -> iterator {
        return {m_cfg, 1u, 0u};
    }

    /// @returns the generator in end state
    DETRAY_HOST_DEVICE
    constexpr auto end() const noexcept -> iterator {
        return {m_cfg, 1u, m_cfg.theta_steps()};
    }

    /// @returns the number of tracks that will be generated
    DETRAY_HOST_DEVICE
    constexpr auto size() const noexcept -> std::size_t {
        return m_cfg.theta_steps() * m_cfg.phi_steps();
    }
};

}  // namespace detray
