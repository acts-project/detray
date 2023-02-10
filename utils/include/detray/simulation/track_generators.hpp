/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/utils/ranges/ranges.hpp"

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
                 point3 trk_origin = {0.f, 0.f, 0.f},
                 scalar trk_mom = 1.f * unit<scalar>::GeV,
                 std::array<scalar, 2> theta_range = {0.01f,
                                                      constant<scalar>::pi},
                 std::array<scalar, 2> phi_range = {-constant<scalar>::pi,
                                                    constant<scalar>::pi},
                 scalar time = 0.f * unit<scalar>::us,
                 scalar charge = -1.f * unit<scalar>::e, std::size_t iph = 1u,
                 std::size_t ith = 0u)
            : m_theta_steps{n_theta},
              m_phi_steps{n_phi},
              m_theta_step_size{(theta_range[1] - theta_range[0]) /
                                static_cast<scalar>(n_theta)},
              m_phi_step_size{(phi_range[1] - phi_range[0]) /
                              static_cast<scalar>(n_phi)},
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
                    m_phi = m_phi_range[0] +
                            static_cast<scalar>(i_phi) * m_phi_step_size;
                    ++i_phi;
                    return *this;
                }
                // Reset phi range
                i_phi = 1;
                m_phi = m_phi_range[0];
                ;
                // Calculate new thetain the given range
                ++i_theta;
                m_theta = m_theta_range[0] +
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
            mom = m_mom_mag * mom;

            return track_t{m_origin, m_time, mom, m_charge};
        }

        /// Start and end values of angle space
        std::size_t m_theta_steps{50u};
        std::size_t m_phi_steps{50u};
        scalar m_theta_step_size{0.f};
        scalar m_phi_step_size{0.f};

        /// Phi and theta angles of momentum direction
        scalar m_phi{-constant<scalar>::pi}, m_theta{0.01f};

        /// Track origin
        point3 m_origin{0.f, 0.f, 0.f};

        /// Magnitude of momentum: Default is one to keep directions normalized
        /// if no momentum information is needed (e.g. for a ray)
        scalar m_mom_mag{1.f * unit<scalar>::GeV};

        /// Range for theta and phi
        std::array<scalar, 2> m_theta_range{0.01f, constant<scalar>::pi};
        std::array<scalar, 2> m_phi_range{-constant<scalar>::pi,
                                          constant<scalar>::pi};

        /// Time parameter and charge of the track
        scalar m_time{0.f}, m_charge{0.f};

        /// Iteration indices
        std::size_t i_phi{0u};
        std::size_t i_theta{0u};
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
    uniform_track_generator(
        std::size_t n_theta, std::size_t n_phi,
        point3 trk_origin = {0.f, 0.f, 0.f},
        scalar trk_mom = 1.f * unit<scalar>::GeV,
        std::array<scalar, 2> theta_range = {0.01f, constant<scalar>::pi},
        std::array<scalar, 2> phi_range = {-constant<scalar>::pi,
                                           constant<scalar>::pi},
        scalar time = 0.f * unit<scalar>::us,
        scalar charge = -1.f * unit<scalar>::e)
        : m_begin{n_theta,   n_phi, trk_origin, trk_mom, theta_range,
                  phi_range, time,  charge,     1u,      0u},
          m_end{n_theta,   n_phi, trk_origin, trk_mom, theta_range,
                phi_range, time,  charge,     1u,      n_theta} {}

    /// Copy assignment operator
    DETRAY_HOST_DEVICE
    uniform_track_generator &operator=(const uniform_track_generator &other) {
        m_begin = other.m_begin;
        m_end = other.m_end;
        return *this;
    }

    /// Move constructor
    uniform_track_generator(uniform_track_generator &&other)
        : m_begin(std::move(other.m_begin)), m_end(std::move(other.m_end)) {}

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
          typename engine_t = std::mt19937_64>
struct random_numbers {

    random_numbers(random_numbers &&other) : engine(std::move(other.engine)) {}

    generator_t gen;
    engine_t engine;

    template <
        typename T = generator_t,
        std::enable_if_t<std::is_same_v<T, std::random_device>, bool> = true>
    random_numbers() : gen{}, engine{gen()} {}

    template <typename T = generator_t,
              std::enable_if_t<std::is_same_v<T, std::seed_seq>, bool> = true>
    random_numbers() : gen{42u}, engine{gen} {}

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

        constexpr iterator() = delete;

        DETRAY_HOST_DEVICE
        iterator(generator_t &rand_gen, std::size_t n_tracks,
                 point3 trk_origin = {0.f, 0.f, 0.f},
                 point3 origin_stddev = {0.f, 0.f, 0.f},
                 std::array<scalar, 2> mom_range = {1.f * unit<scalar>::GeV,
                                                    1.f * unit<scalar>::GeV},
                 std::array<scalar, 2> theta_range = {0.01f,
                                                      constant<scalar>::pi},
                 std::array<scalar, 2> phi_range = {-constant<scalar>::pi,
                                                    constant<scalar>::pi},
                 scalar time = 0.f * unit<scalar>::us,
                 scalar charge = -1.f * unit<scalar>::e)
            : m_rnd_numbers{rand_gen},
              m_tracks{n_tracks},
              m_origin{trk_origin},
              m_origin_stddev{origin_stddev},
              m_mom_mag{rand_gen(mom_range[0], mom_range[1])},
              m_phi{rand_gen(phi_range[0], phi_range[1])},
              m_theta{rand_gen(theta_range[0], theta_range[1])},
              m_mom_range{mom_range},
              m_phi_range{phi_range},
              m_theta_range{theta_range},
              m_time{time},
              m_charge{charge} {
            m_vertex = {
                std::normal_distribution<scalar>(
                    trk_origin[0], origin_stddev[0])(m_rnd_numbers.engine),
                std::normal_distribution<scalar>(
                    trk_origin[1], origin_stddev[1])(m_rnd_numbers.engine),
                std::normal_distribution<scalar>(
                    trk_origin[2], origin_stddev[2])(m_rnd_numbers.engine)};
        }

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
            m_vertex = {
                std::normal_distribution<scalar>(
                    m_origin[0], m_origin_stddev[0])(m_rnd_numbers.engine),
                std::normal_distribution<scalar>(
                    m_origin[1], m_origin_stddev[1])(m_rnd_numbers.engine),
                std::normal_distribution<scalar>(
                    m_origin[2], m_origin_stddev[2])(m_rnd_numbers.engine)};
            m_mom_mag = std::max(m_rnd_numbers(m_mom_range[0], m_mom_range[1]),
                                 scalar(0.f));
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
                        std::sin(m_phi) * std::sin(m_theta),
                        math_ns::cos(m_theta)};
            // Magnitude of momentum
            vector::normalize(mom);
            mom = m_mom_mag * mom;

            return track_t{m_vertex, m_time, mom, m_charge};
        }

        /// Random number generator
        generator_t &m_rnd_numbers;

        /// How many tracks will be generated
        std::size_t m_tracks{0u};

        /// Track origin
        point3 m_origin{0.f, 0.f, 0.f};

        /// Track origin standard deviatation
        point3 m_origin_stddev{0.f, 0.f, 0.f};

        /// Track vertex
        point3 m_vertex{0.f, 0.f, 0.f};

        /// Magnitude of momentum: Default is one to keep directions normalized
        /// if no momentum information is needed (e.g. for a ray)
        scalar m_mom_mag{1.f * unit<scalar>::GeV};

        /// Phi and theta angles of momentum direction (random)
        scalar m_phi{-constant<scalar>::pi}, m_theta{0.01f};

        /// Range for theta and phi
        std::array<scalar, 2> m_mom_range{1.f * unit<scalar>::GeV,
                                          1.f * unit<scalar>::GeV};
        std::array<scalar, 2> m_phi_range{-constant<scalar>::pi,
                                          constant<scalar>::pi};
        std::array<scalar, 2> m_theta_range{0.01f, constant<scalar>::pi};

        /// Time parameter and charge of the track
        scalar m_time{0.f}, m_charge{0.f};
    };

    generator_t m_gen;
    iterator m_begin{}, m_end{};

    public:
    using iterator_t = iterator;

    /// Default constructor
    constexpr random_track_generator() = default;

    /// Move constructor
    random_track_generator(random_track_generator &&other)
        : m_gen(std::move(other.m_gen)),
          m_begin(std::move(other.m_begin)),
          m_end(std::move(other.m_end)) {}

    /// Paramtetrized constructor for fine-grained configurations
    ///
    /// @param n_tracks the number of steps in the theta space
    /// @param trk_origin the starting point of the track
    /// @param origin_stddev the standard deviation of origin
    /// @param mom_range the range of the track momentum (in GeV)
    /// @param theta_range the range for theta values
    /// @param phi_range the range for phi values
    /// @param time time measurement (micro seconds)
    /// @param charge charge of particle (e)
    DETRAY_HOST_DEVICE
    random_track_generator(
        std::size_t n_tracks, point3 trk_origin = {0.f, 0.f, 0.f},
        point3 origin_stddev = {0.f, 0.f, 0.f},
        std::array<scalar, 2> mom_range = {1.f * unit<scalar>::GeV,
                                           1.f * unit<scalar>::GeV},
        std::array<scalar, 2> theta_range = {0.01f, constant<scalar>::pi},
        std::array<scalar, 2> phi_range = {-constant<scalar>::pi,
                                           constant<scalar>::pi},
        scalar time = 0.f * unit<scalar>::us,
        scalar charge = -1.f * unit<scalar>::e)
        : m_gen{},
          m_begin{m_gen,       0u,        trk_origin, origin_stddev, mom_range,
                  theta_range, phi_range, time,       charge},
          m_end{m_gen,       n_tracks,  trk_origin, origin_stddev, mom_range,
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
