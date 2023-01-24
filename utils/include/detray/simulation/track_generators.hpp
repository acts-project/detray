/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/utils/ranges/ranges.hpp"
#include "detray/definitions/math.hpp"

// System include(s)
#include <algorithm>
#include <cmath>
#include <random>

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
    private:
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /// @brief Nested iterator type that generates track states.
    struct iterator {

        using difference_type = std::ptrdiff_t;
        using value_type = track_t;
        using pointer = track_t *;
        using reference = track_t &;
        using iterator_category = detray::ranges::input_iterator_tag;

        constexpr iterator() = default;

        DETRAY_HOST_DEVICE
        iterator(std::size_t n_theta, std::size_t n_phi,
                 point3 trk_origin = {0., 0., 0.},
                 scalar trk_mom = 1. * unit<scalar>::GeV,
                 std::array<scalar, 2> theta_range = {0.01, M_PI},
                 std::array<scalar, 2> phi_range = {-M_PI, M_PI},
                 scalar time = 0. * unit<scalar>::us,
                 scalar charge = -1. * unit<scalar>::e, std::size_t iph = 1,
                 std::size_t ith = 0)
            : m_theta_steps{n_theta},
              m_phi_steps{n_phi},
              m_theta_step_size{(theta_range[1] - theta_range[0]) / n_theta},
              m_phi_step_size{(phi_range[1] - phi_range[0]) / n_phi},
              m_phi{phi_range[0]},
              m_theta{theta_range[0]},
              m_origin{trk_origin},
              m_mom_mag{trk_mom},
              m_theta_range{theta_range},
              m_phi_range{phi_range},
              m_time{time},
              m_charge{charge},
              i_phi{iph},
              i_theta{ith} {}

        /// @returns whether we reached end of angle space
        DETRAY_HOST_DEVICE
        constexpr bool operator==(const iterator &rhs) const {
            return rhs.m_theta_steps == m_theta_steps and
                   rhs.m_phi_steps == m_phi_steps and rhs.i_phi == i_phi and
                   rhs.i_theta == i_theta;
        }

        /// @returns whether we reached end of angle space
        DETRAY_HOST_DEVICE
        constexpr bool operator!=(const iterator &rhs) const {
            return not(*this == rhs);
        }

        /// Iterate through angle space according to given step sizes.
        ///
        /// @returns the generator at its next position.
        DETRAY_HOST_DEVICE
        constexpr auto operator++() -> iterator & {
            // Check theta range according to step size
            if (i_theta < m_theta_steps) {
                // Check phi sub-range
                if (i_phi < m_phi_steps) {
                    // Calculate new phi in the given range
                    m_phi = m_phi_range[0] + i_phi * m_phi_step_size;
                    ++i_phi;
                    return *this;
                }
                // Reset phi range
                i_phi = 1;
                m_phi = m_phi_range[0];
                ;
                // Calculate new thetain the given range
                ++i_theta;
                m_theta = m_theta_range[0] + i_theta * m_theta_step_size;
            }
            return *this;
        }

        /// @returns a track instance from generated momentum direction
        DETRAY_HOST_DEVICE
        track_t operator*() const {
            // Momentum direction from angles
            vector3 mom{math_ns::cos(m_phi) * std::sin(m_theta),
                        std::sin(m_phi) * std::sin(m_theta), math_ns::cos(m_theta)};
            // Magnitude of momentum
            vector::normalize(mom);
            mom = m_mom_mag * mom;

            return track_t{m_origin, m_time, mom, m_charge};
        }

        /// Start and end values of angle space
        std::size_t m_theta_steps{50};
        std::size_t m_phi_steps{50};
        scalar m_theta_step_size{0};
        scalar m_phi_step_size{0};

        /// Phi and theta angles of momentum direction
        scalar m_phi{-M_PI}, m_theta{0.01};

        /// Track origin
        point3 m_origin{0., 0., 0.};

        /// Magnitude of momentum: Default is one to keep directions normalized
        /// if no momentum information is needed (e.g. for a ray)
        scalar m_mom_mag{1. * unit<scalar>::GeV};

        /// Range for theta and phi
        std::array<scalar, 2> m_theta_range{0.01, M_PI};
        std::array<scalar, 2> m_phi_range{-M_PI, M_PI};

        /// Time parameter and charge of the track
        scalar m_time{0}, m_charge{0};

        /// Iteration indices
        std::size_t i_phi{0};
        std::size_t i_theta{0};
    };

    iterator m_begin{}, m_end{};

    public:
    using iterator_t = iterator;

    /// Default constructor
    constexpr uniform_track_generator() = default;

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
    uniform_track_generator(std::size_t n_theta, std::size_t n_phi,
                            point3 trk_origin = {0., 0., 0.},
                            scalar trk_mom = 1. * unit<scalar>::GeV,
                            std::array<scalar, 2> theta_range = {0.01, M_PI},
                            std::array<scalar, 2> phi_range = {-M_PI, M_PI},
                            scalar time = 0. * unit<scalar>::us,
                            scalar charge = -1. * unit<scalar>::e)
        : m_begin{n_theta,   n_phi, trk_origin, trk_mom, theta_range,
                  phi_range, time,  charge,     1,       0},
          m_end{n_theta,   n_phi, trk_origin, trk_mom, theta_range,
                phi_range, time,  charge,     1,       n_theta} {}

    /// Copy assignment operator
    DETRAY_HOST_DEVICE
    uniform_track_generator &operator=(const uniform_track_generator &other) {
        m_begin = other.m_begin;
        m_end = other.m_end;
        return *this;
    }

    /// @returns the generator in initial state: Default values reflect the
    /// first phi angle iteration.
    DETRAY_HOST_DEVICE
    constexpr auto begin() const noexcept -> iterator { return m_begin; }

    /// @returns the generator in end state
    DETRAY_HOST_DEVICE
    constexpr auto end() const noexcept -> iterator { return m_end; }

    /// @returns the number of tracks that will be generated
    DETRAY_HOST_DEVICE
    constexpr auto size() const noexcept -> std::size_t {
        return m_begin.m_theta_steps * m_begin.m_phi_steps;
    }
};

/// Wrapper for random number generatrion for the @c random_track_generator
template <typename scalar_t = scalar,
          typename distribution_t = std::uniform_real_distribution<scalar_t>,
          typename generator_t = std::random_device,
          typename engine_t = std::mt19937>
struct random_numbers {
    generator_t gen;
    engine_t engine;

    template <
        typename T = generator_t,
        std::enable_if_t<std::is_same_v<T, std::random_device>, bool> = true>
    random_numbers() : gen{}, engine{gen()} {}

    template <typename T = generator_t,
              std::enable_if_t<std::is_same_v<T, std::seed_seq>, bool> = true>
    random_numbers() : gen{42}, engine{gen} {}

    template <typename T = distribution_t,
              std::enable_if_t<
                  std::is_same_v<T, std::uniform_real_distribution<scalar_t>>,
                  bool> = true>
    DETRAY_HOST auto operator()(const scalar_t min, const scalar_t max) {
        return distribution_t(min, max)(engine);
    }

    template <
        typename T = distribution_t,
        std::enable_if_t<std::is_same_v<T, std::normal_distribution<scalar_t>>,
                         bool> = true>
    DETRAY_HOST auto operator()(const scalar_t min, const scalar_t max) {
        scalar_t mu{min + 0.5f * (max - min)};
        return distribution_t(mu, 0.5f / 3.0f * (max - min))(engine);
    }
};

/// @brief Generates track states with random momentum directions.
///
/// Generates the phi and theta angles of the track momentum according to a
/// given random number distribution.
///
/// @tparam track_t the type of track parametrization that should be used.
/// @tparam generator_t source of random numbers
///
/// @note Since the random number generator might not be copy constructible,
/// neither is this generator. The iterators hold a reference to the rand
/// generator, which must not be invalidated during the iteration.
/// @note the random numbers are clamped to fit the phi/theta ranges. This can
/// effect distribution mean etc.
template <typename track_t, typename generator_t = random_numbers<>>
class random_track_generator
    : public detray::ranges::view_interface<random_track_generator<track_t>> {
    private:
    using point3 = __plugin::point3<detray::scalar>;
    using vector3 = __plugin::vector3<detray::scalar>;

    /// @brief Nested iterator type that generates track states.
    struct iterator {

        using difference_type = std::ptrdiff_t;
        using value_type = track_t;
        using pointer = track_t *;
        using reference = track_t &;
        using iterator_category = detray::ranges::input_iterator_tag;

        iterator() = delete;

        DETRAY_HOST_DEVICE
        iterator(generator_t &rand_gen, std::size_t n_tracks,
                 point3 trk_origin = {0., 0., 0.},
                 scalar trk_mom = 1. * unit<scalar>::GeV,
                 std::array<scalar, 2> theta_range = {0.01, M_PI},
                 std::array<scalar, 2> phi_range = {-M_PI, M_PI},
                 scalar time = 0. * unit<scalar>::us,
                 scalar charge = -1. * unit<scalar>::e)
            : m_rnd_numbers{rand_gen},
              m_tracks{n_tracks},
              m_origin{trk_origin},
              m_mom_mag{trk_mom},
              m_phi{rand_gen(phi_range[0], phi_range[1])},
              m_theta{rand_gen(theta_range[0], theta_range[1])},
              m_phi_range{phi_range},
              m_theta_range{theta_range},
              m_time{time},
              m_charge{charge} {}

        /// @returns whether we reached the end of iteration
        DETRAY_HOST_DEVICE
        constexpr bool operator==(const iterator &rhs) const {
            return rhs.m_tracks == m_tracks;
        }

        /// @returns whether we reached the end of iteration
        DETRAY_HOST_DEVICE
        constexpr bool operator!=(const iterator &rhs) const {
            return not(*this == rhs);
        }

        /// Genrate a new random direction.
        ///
        /// @returns the generator at its next position.
        DETRAY_HOST_DEVICE
        auto operator++() -> iterator & {
            m_phi = std::clamp(m_rnd_numbers(m_phi_range[0], m_phi_range[1]),
                               m_phi_range[0], m_phi_range[1]);
            m_theta =
                std::clamp(m_rnd_numbers(m_theta_range[0], m_theta_range[1]),
                           m_theta_range[0], m_theta_range[1]);

            ++m_tracks;

            return *this;
        }

        /// @returns a track instance from generated momentum direction
        DETRAY_HOST_DEVICE
        track_t operator*() const {
            // Momentum direction from angles
            vector3 mom{math_ns::cos(m_phi) * std::sin(m_theta),
                        std::sin(m_phi) * std::sin(m_theta), math_ns::cos(m_theta)};
            // Magnitude of momentum
            vector::normalize(mom);
            mom = m_mom_mag * mom;

            return track_t{m_origin, m_time, mom, m_charge};
        }

        /// Random number generator
        generator_t &m_rnd_numbers;

        /// How many tracks will be generated
        std::size_t m_tracks{0};

        /// Track origin
        point3 m_origin{0., 0., 0.};

        /// Magnitude of momentum: Default is one to keep directions normalized
        /// if no momentum information is needed (e.g. for a ray)
        scalar m_mom_mag{1. * unit<scalar>::GeV};

        /// Phi and theta angles of momentum direction (random)
        scalar m_phi{-M_PI}, m_theta{0.01};

        /// Range for theta and phi
        std::array<scalar, 2> m_phi_range{-M_PI, M_PI};
        std::array<scalar, 2> m_theta_range{0.01, M_PI};

        /// Time parameter and charge of the track
        scalar m_time{0}, m_charge{0};
    };

    generator_t m_gen;
    iterator m_begin{}, m_end{};

    public:
    using iterator_t = iterator;

    /// Default constructor
    random_track_generator() = default;

    /// Paramtetrized constructor for fine-grained configurations
    ///
    /// @param n_tracks the number of steps in the theta space
    /// @param trk_origin the starting point of the track
    /// @param trk_mom magnitude of the track momentum (in GeV)
    /// @param theta_range the range for theta values
    /// @param phi_range the range for phi values
    /// @param time time measurement (micro seconds)
    /// @param charge charge of particle (e)
    DETRAY_HOST_DEVICE
    random_track_generator(std::size_t n_tracks,
                           point3 trk_origin = {0., 0., 0.},
                           scalar trk_mom = 1. * unit<scalar>::GeV,
                           std::array<scalar, 2> theta_range = {0.01, M_PI},
                           std::array<scalar, 2> phi_range = {-M_PI, M_PI},
                           scalar time = 0. * unit<scalar>::us,
                           scalar charge = -1. * unit<scalar>::e)
        : m_gen{},
          m_begin{m_gen,       0,         trk_origin, trk_mom,
                  theta_range, phi_range, time,       charge},
          m_end{m_gen,       n_tracks,  trk_origin, trk_mom,
                theta_range, phi_range, time,       charge} {}

    /// @returns the generator in initial state.
    /// @note the underlying random number generator has deleted copy
    /// constructor, so the iterator needs to be built from scratch
    DETRAY_HOST_DEVICE
    auto begin() const noexcept -> iterator {
        // return {0, m_begin.m_origin, m_begin.m_mom_mag, m_begin.m_phi_range,
        // m_begin.m_theta_range, m_begin.m_time, m_begin.m_charge};
        return m_begin;
    }

    /// @returns the generator in end state
    DETRAY_HOST_DEVICE
    constexpr auto end() const noexcept -> iterator { return m_end; }

    /// @returns the number of tracks that will be generated
    DETRAY_HOST_DEVICE
    constexpr auto size() const noexcept -> std::size_t {
        return m_end.m_tracks;
    }
};

}  // namespace detray
