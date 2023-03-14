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

/// Wrapper for random number generatrion for the @c random_track_generator
template <typename scalar_t = scalar,
          typename distribution_t = std::uniform_real_distribution<scalar_t>,
          typename generator_t = std::random_device,
          typename engine_t = std::mt19937_64>
struct random_numbers {

    random_numbers(random_numbers&& other) : engine(std::move(other.engine)) {}

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
    public:
    using point3 = typename track_t::point3;
    using vector3 = typename track_t::vector3;

    /// Configure how tracks are generated
    struct configuration {

        /// Gaussian vertex smearing
        bool m_do_vtx_smearing = true;

        /// How many tracks will be generated
        std::size_t m_n_tracks{10u};

        /// Range for theta and phi
        std::array<scalar, 2> m_phi_range{-constant<scalar>::pi,
                                          constant<scalar>::pi};
        std::array<scalar, 2> m_theta_range{0.01f, constant<scalar>::pi};
        /// Momentum range
        std::array<scalar, 2> m_mom_range{1.f * unit<scalar>::GeV,
                                          1.f * unit<scalar>::GeV};
        /// Track origin
        point3 m_origin{0.f, 0.f, 0.f}, m_origin_stddev{0.f, 0.f, 0.f};

        /// Time parameter and charge of the track
        scalar m_time{0.f * unit<scalar>::us}, m_charge{-1.f * unit<scalar>::e};

        /// Setters
        /// @{
        configuration& do_vertex_smearing(bool b) {
            m_do_vtx_smearing = b;
            return *this;
        }
        configuration& n_tracks(std::size_t n) {
            m_n_tracks = n;
            return *this;
        }
        configuration& theta_range(scalar low, scalar high) {
            m_theta_range = {low, high};
            return *this;
        }
        configuration& phi_range(scalar low, scalar high) {
            m_phi_range = {low, high};
            return *this;
        }
        configuration& origin(point3 ori) {
            m_origin = ori;
            return *this;
        }
        configuration& origin_stddev(point3 stddev) {
            m_origin_stddev = stddev;
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
        constexpr bool do_vertex_smearing() const { return m_do_vtx_smearing; }
        constexpr std::size_t n_tracks() const { return m_n_tracks; }
        constexpr const std::array<scalar, 2>& theta_range() const {
            return m_theta_range;
        }
        constexpr const std::array<scalar, 2>& phi_range() const {
            return m_phi_range;
        }
        constexpr const std::array<scalar, 2>& mom_range() const {
            return m_mom_range;
        }
        constexpr const point3& origin() const { return m_origin; }
        constexpr const point3& origin_stddev() const {
            return m_origin_stddev;
        }
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

        constexpr iterator() = delete;

        DETRAY_HOST_DEVICE
        iterator(generator_t& rand_gen, configuration cfg, std::size_t n_tracks)
            : m_rnd_numbers{rand_gen}, m_tracks{n_tracks}, m_cfg{cfg} {}

        /// @returns whether we reached the end of iteration
        DETRAY_HOST_DEVICE
        constexpr bool operator==(const iterator& rhs) const {
            return rhs.m_tracks == m_tracks;
        }

        /// @returns whether we reached the end of iteration
        DETRAY_HOST_DEVICE
        constexpr bool operator!=(const iterator& rhs) const {
            return not(*this == rhs);
        }

        /// @returns the generator at its next position.
        DETRAY_HOST_DEVICE
        constexpr auto operator++() -> iterator& {
            ++m_tracks;
            return *this;
        }

        /// @returns a track instance from random-generated momentum
        DETRAY_HOST_DEVICE
        track_t operator*() const {

            const point3 vtx =
                m_cfg.do_vertex_smearing()
                    ? point3{std::normal_distribution<scalar>(
                                 m_cfg.origin()[0], m_cfg.origin_stddev()[0])(
                                 m_rnd_numbers.engine),
                             std::normal_distribution<scalar>(
                                 m_cfg.origin()[1], m_cfg.origin_stddev()[1])(
                                 m_rnd_numbers.engine),
                             std::normal_distribution<scalar>(
                                 m_cfg.origin()[2], m_cfg.origin_stddev()[2])(
                                 m_rnd_numbers.engine)}
                    : m_cfg.origin();

            const std::array<scalar, 2>& phi_rng = m_cfg.phi_range();
            const std::array<scalar, 2>& theta_rng = m_cfg.theta_range();
            const std::array<scalar, 2>& mom_rng = m_cfg.mom_range();
            scalar p_mag{
                std::max(m_rnd_numbers(mom_rng[0], mom_rng[1]), scalar{0.f})};
            scalar phi{std::clamp(m_rnd_numbers(phi_rng[0], phi_rng[1]),
                                  phi_rng[0], phi_rng[1])};
            scalar theta{std::clamp(m_rnd_numbers(theta_rng[0], theta_rng[1]),
                                    theta_rng[0], theta_rng[1])};

            // Momentum direction from angles
            vector3 mom{math_ns::cos(phi) * math_ns::sin(theta),
                        math_ns::sin(phi) * math_ns::sin(theta),
                        math_ns::cos(theta)};
            // Magnitude of momentum
            vector::normalize(mom);
            mom = p_mag * mom;

            return track_t{vtx, m_cfg.time(), mom, m_cfg.charge()};
        }

        /// Random number generator
        generator_t& m_rnd_numbers;

        /// How many tracks will be generated
        std::size_t m_tracks{0u};

        /// Configuration
        configuration m_cfg{};
    };

    generator_t m_gen;
    configuration m_cfg{};

    public:
    using iterator_t = iterator;

    /// Default constructor
    constexpr random_track_generator() = default;

    /// Construct from external configuration
    constexpr random_track_generator(configuration cfg) : m_gen{}, m_cfg(cfg) {}

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
        scalar charge = -1.f * unit<scalar>::e, bool do_vertex_smearing = true)
        : m_gen{}, m_cfg{do_vertex_smearing, n_tracks,  phi_range,
                         theta_range,        mom_range, trk_origin,
                         origin_stddev,      time,      charge} {}

    /// Move constructor
    random_track_generator(random_track_generator&& other)
        : m_gen(std::move(other.m_gen)), m_cfg(std::move(other.m_cfg)) {}

    /// Access the configuration
    constexpr configuration& config() { return m_cfg; }

    /// @returns the generator in initial state.
    /// @note the underlying random number generator has deleted copy
    /// constructor, so the iterator needs to be built from scratch
    DETRAY_HOST_DEVICE
    auto begin() noexcept -> iterator { return {m_gen, m_cfg, 0u}; }

    /// @returns the generator in end state
    DETRAY_HOST_DEVICE
    auto end() noexcept -> iterator { return {m_gen, m_cfg, m_cfg.n_tracks()}; }

    /// @returns the number of tracks that will be generated
    DETRAY_HOST_DEVICE
    constexpr auto size() const noexcept -> std::size_t {
        return m_cfg.n_tracks();
    }
};

}  // namespace detray