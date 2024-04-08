/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
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
template <typename track_t>
class uniform_track_generator
    : public detray::ranges::view_interface<uniform_track_generator<track_t>> {

    using point3 = typename track_t::point3_type;
    using vector3 = typename track_t::vector3_type;

    public:
    using track_type = track_t;

    /// Configure how tracks are generated
    struct configuration {
        /// Ensure sensible values at the theta bounds, even in single precision
        static constexpr scalar epsilon{1e-2f};

        /// Range for theta and phi
        std::array<scalar, 2> m_phi_range{-constant<scalar>::pi,
                                          constant<scalar>::pi};
        std::array<scalar, 2> m_theta_range{epsilon,
                                            constant<scalar>::pi - epsilon};
        std::array<scalar, 2> m_eta_range{-5.f, 5.f};

        /// Angular step size
        std::size_t m_phi_steps{50u};
        std::size_t m_theta_steps{50u};

        /// Do uniform eta steps instead of uniform theta steps
        /// (use same number of steps and range)
        bool m_uniform_eta{false};

        /// Track origin
        point3 m_origin{0.f, 0.f, 0.f};

        /// Magnitude of momentum: Default is one to keep directions normalized
        /// if no momentum information is needed (e.g. for a ray)
        scalar m_p_mag{1.f * unit<scalar>::GeV};
        /// Whether to interpret the momentum @c m_p_mag as p_T
        bool m_is_pT{false};

        /// Time parameter and charge of the track
        scalar m_time{0.f * unit<scalar>::us};
        scalar m_charge{-1.f * unit<scalar>::e};

        /// Setters
        /// @{
        DETRAY_HOST_DEVICE configuration& phi_range(scalar low, scalar high) {
            assert(low <= high);
            m_phi_range = {low, high};
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& theta_range(scalar low, scalar high) {
            auto min_theta{
                std::clamp(low, epsilon, constant<scalar>::pi - epsilon)};
            auto max_theta{
                std::clamp(high, epsilon, constant<scalar>::pi - epsilon)};

            assert(min_theta <= max_theta);

            m_theta_range = {min_theta, max_theta};
            m_uniform_eta = false;
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& eta_range(scalar low, scalar high) {
            // This value is more or less random
            constexpr auto num_max{0.001f * std::numeric_limits<scalar>::max()};
            auto min_eta{low > -num_max ? low : -num_max};
            auto max_eta{high < num_max ? high : num_max};

            assert(min_eta <= max_eta);

            m_eta_range = {min_eta, max_eta};
            m_uniform_eta = true;
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& phi_steps(std::size_t n) {
            assert(n > 0);
            m_phi_steps = n;
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& theta_steps(std::size_t n) {
            assert(n > 0);
            m_theta_steps = n;
            m_uniform_eta = false;
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& eta_steps(std::size_t n) {
            assert(n > 0);
            m_theta_steps = n;
            m_uniform_eta = true;
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& uniform_eta(bool b) {
            m_uniform_eta = b;
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& origin(point3 ori) {
            m_origin = ori;
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& p_tot(scalar p) {
            m_is_pT = false;
            m_p_mag = p;
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& p_T(scalar p) {
            m_is_pT = true;
            m_p_mag = p;
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& time(scalar t) {
            m_time = t;
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& charge(scalar q) {
            m_charge = q;
            return *this;
        }
        /// @}

        /// Getters
        /// @{
        DETRAY_HOST_DEVICE constexpr std::array<scalar, 2> phi_range() const {
            return m_phi_range;
        }
        DETRAY_HOST_DEVICE constexpr std::array<scalar, 2> theta_range() const {
            return m_theta_range;
        }
        DETRAY_HOST_DEVICE constexpr std::array<scalar, 2> eta_range() const {
            return m_eta_range;
        }
        DETRAY_HOST_DEVICE constexpr std::size_t phi_steps() const {
            return m_phi_steps;
        }
        DETRAY_HOST_DEVICE constexpr std::size_t theta_steps() const {
            return m_theta_steps;
        }
        DETRAY_HOST_DEVICE constexpr std::size_t eta_steps() const {
            return m_theta_steps;
        }
        DETRAY_HOST_DEVICE constexpr bool uniform_eta() const {
            return m_uniform_eta;
        }
        DETRAY_HOST_DEVICE constexpr const point3& origin() const {
            return m_origin;
        }
        DETRAY_HOST_DEVICE constexpr bool is_pT() const { return m_is_pT; }
        DETRAY_HOST_DEVICE constexpr scalar time() const { return m_time; }
        DETRAY_HOST_DEVICE constexpr scalar charge() const { return m_charge; }
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
              m_phi_step_size{(cfg.phi_range()[1] - cfg.phi_range()[0]) /
                              static_cast<scalar>(cfg.phi_steps())},
              m_theta_step_size{(cfg.theta_range()[1] - cfg.theta_range()[0]) /
                                static_cast<scalar>(cfg.theta_steps() - 1u)},
              m_eta_step_size{(cfg.eta_range()[1] - cfg.eta_range()[0]) /
                              static_cast<scalar>(cfg.eta_steps() - 1u)},
              m_phi{cfg.phi_range()[0]},
              m_theta{cfg.uniform_eta() ? get_theta(cfg.eta_range()[0])
                                        : cfg.theta_range()[0]},
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

                if (m_cfg.uniform_eta()) {
                    const scalar eta =
                        m_cfg.eta_range()[0] +
                        static_cast<scalar>(i_theta) * m_eta_step_size;
                    m_theta = get_theta(eta);
                } else {
                    m_theta = m_cfg.theta_range()[0] +
                              static_cast<scalar>(i_theta) * m_theta_step_size;
                }
            }
            return *this;
        }

        /// @returns a track instance from generated momentum direction
        DETRAY_HOST_DEVICE
        track_t operator*() const {

            const scalar sin_theta{math::sin(m_theta)};

            // Momentum direction from angles
            vector3 p{math::cos(m_phi) * sin_theta,
                      math::sin(m_phi) * sin_theta, math::cos(m_theta)};
            // Magnitude of momentum
            vector::normalize(p);

            p = (m_cfg.is_pT() ? 1.f / sin_theta : 1.f) * m_cfg.m_p_mag * p;

            return track_t{m_cfg.origin(), m_cfg.time(), p, m_cfg.charge()};
        }

        /// Current configuration
        configuration m_cfg{};

        /// Angular step sizes
        scalar m_phi_step_size{0.f};
        scalar m_theta_step_size{0.f};
        scalar m_eta_step_size{0.f};

        /// Phi and theta angles of momentum direction
        scalar m_phi{-constant<scalar>::pi};
        scalar m_theta{configuration::epsilon};

        /// Iteration indices
        std::size_t i_phi{0u};
        std::size_t i_theta{0u};

        private:
        /// @returns the theta angle for a given @param eta value
        DETRAY_HOST_DEVICE
        scalar get_theta(const scalar eta) {
            return 2.f * math::atan(math::exp(-eta));
        }
    };

    configuration m_cfg{};

    public:
    using iterator_t = iterator;

    /// Default constructor
    constexpr uniform_track_generator() = default;

    /// Construct from external configuration @param cfg
    DETRAY_HOST_DEVICE
    constexpr uniform_track_generator(configuration cfg) : m_cfg{cfg} {}

    /// Paramtetrized constructor for quick construction of simple tasks
    ///
    /// @note For more complex tasks, use the @c configuration type
    ///
    /// @param n_theta the number of steps in the theta space
    /// @param n_phi the number of steps in the phi space
    /// @param p_mag magnitude of the track momentum (in GeV)
    /// @param uniform_eta uniformly step through eta space instead of theta
    /// @param charge charge of particle (e)
    DETRAY_HOST_DEVICE
    uniform_track_generator(std::size_t n_phi, std::size_t n_theta,
                            scalar p_mag = 1.f * unit<scalar>::GeV,
                            bool uniform_eta = false,
                            scalar charge = -1.f * unit<scalar>::e)
        : m_cfg{} {
        m_cfg.phi_steps(n_phi).theta_steps(n_theta);
        m_cfg.uniform_eta(uniform_eta);
        m_cfg.p_tot(p_mag);
        m_cfg.charge(charge);
    }

    /// Move constructor
    DETRAY_HOST_DEVICE
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
    DETRAY_HOST_DEVICE
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
        return m_cfg.phi_steps() * m_cfg.theta_steps();
    }
};

}  // namespace detray
