/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/utils/ranges/ranges.hpp"

// System include(s)
#include <algorithm>
#include <array>
#include <random>

namespace detray {

/// Wrapper for CPU random number generatrion for the @c random_track_generator
template <typename scalar_t = scalar,
          typename distribution_t = std::uniform_real_distribution<scalar_t>,
          typename engine_t = std::mt19937_64>
struct random_numbers {

    std::seed_seq m_seeds;
    engine_t m_engine;

    /// Random seed for increased entropy. Different random numbers each time
    DETRAY_HOST
    random_numbers() : m_seeds{get_random_seeds()}, m_engine{m_seeds} {}

    /// Deterministic @param seed for reproducible results
    DETRAY_HOST
    random_numbers(std::size_t seed) : m_seeds{seed}, m_engine{m_seeds} {}

    /// Copy constructor
    DETRAY_HOST
    random_numbers(random_numbers&& other)
        : m_engine(std::move(other.m_engine)) {}

    /// Generate random numbers in a given range
    DETRAY_HOST auto operator()(const std::array<scalar_t, 2> range = {
                                    -std::numeric_limits<scalar_t>::max(),
                                    std::numeric_limits<scalar_t>::max()}) {
        const scalar_t min{range[0]}, max{range[1]};
        assert(min <= max);

        // Uniform
        if constexpr (std::is_same_v<
                          distribution_t,
                          std::uniform_real_distribution<scalar_t>>) {

            return distribution_t(min, max)(m_engine);

            // Normal
        } else if constexpr (std::is_same_v<
                                 distribution_t,
                                 std::normal_distribution<scalar_t>>) {

            scalar_t mu{min + 0.5f * (max - min)};
            return distribution_t(mu, 0.5f / 3.0f * (max - min))(m_engine);
        }
    }

    /// Explicit normal distribution around a @param mean and @param stddev
    DETRAY_HOST auto normal(const scalar_t mean, const scalar_t stddev) {
        return std::normal_distribution<scalar_t>(mean, stddev)(m_engine);
    }

    /// Get proper random seeds for RNG
    DETRAY_HOST
    auto get_random_seeds() {
        std::random_device rd{};
        if constexpr (std::is_same_v<engine_t, std::mt19937> ||
                      std::is_same_v<engine_t, std::mt19937_64>) {
            // Mersenne Twister 19937 needs 624(312) 32(64)-bit seeds to fully
            // initialize its state
            std::array<typename engine_t::result_type, engine_t::state_size>
                seeds;
            std::generate_n(seeds.begin(), seeds.size(), std::ref(rd));
            return std::seed_seq{seeds.begin(), seeds.end()};
        } else {
            return std::seed_seq{rd()};
        }
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
    : public detray::ranges::view_interface<
          random_track_generator<track_t, generator_t>> {
    public:
    using point3 = typename track_t::point3;
    using vector3 = typename track_t::vector3;

    /// Configure how tracks are generated
    struct configuration {
        /// Ensure sensible values at the theta bounds, even in single precision
        static constexpr scalar epsilon{1e-2f};

        /// Gaussian vertex smearing
        bool m_do_vtx_smearing = true;

        /// Monte-Carlo seed
        std::size_t m_seed{0u};

        /// How many tracks will be generated
        std::size_t m_n_tracks{10u};

        /// Range for phi and theta
        std::array<scalar, 2> m_phi_range{-constant<scalar>::pi,
                                          constant<scalar>::pi};
        std::array<scalar, 2> m_theta_range{epsilon,
                                            constant<scalar>::pi - epsilon};

        /// Momentum range
        std::array<scalar, 2> m_mom_range{1.f * unit<scalar>::GeV,
                                          1.f * unit<scalar>::GeV};
        /// Whether to interpret the momentum @c m_mom_range as p_T
        bool m_is_pT{false};

        /// Track origin
        point3 m_origin{0.f, 0.f, 0.f}, m_origin_stddev{0.f, 0.f, 0.f};

        /// Time parameter and charge of the track
        scalar m_time{0.f * unit<scalar>::us}, m_charge{-1.f * unit<scalar>::e};

        /// Setters
        /// @{
        DETRAY_HOST_DEVICE configuration& seed(const std::size_t s) {
            m_seed = s;
            return *this;
        }

        DETRAY_HOST_DEVICE configuration& do_vertex_smearing(bool b) {
            m_do_vtx_smearing = b;
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& n_tracks(std::size_t n) {
            assert(n > 0);
            m_n_tracks = n;
            return *this;
        }
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
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& mom_range(scalar low, scalar high) {
            m_is_pT = false;
            assert(low <= high);
            m_mom_range = {low, high};
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& pT_range(scalar low, scalar high) {
            m_is_pT = true;
            assert(low <= high);
            m_mom_range = {low, high};
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& origin(point3 ori) {
            m_origin = ori;
            return *this;
        }
        DETRAY_HOST_DEVICE configuration& origin_stddev(point3 stddev) {
            m_origin_stddev = stddev;
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
        DETRAY_HOST_DEVICE constexpr std::size_t seed() const { return m_seed; }
        DETRAY_HOST_DEVICE constexpr bool do_vertex_smearing() const {
            return m_do_vtx_smearing;
        }
        DETRAY_HOST_DEVICE constexpr std::size_t n_tracks() const {
            return m_n_tracks;
        }
        DETRAY_HOST_DEVICE constexpr const std::array<scalar, 2>& phi_range()
            const {
            return m_phi_range;
        }
        DETRAY_HOST_DEVICE constexpr const std::array<scalar, 2>& theta_range()
            const {
            return m_theta_range;
        }
        DETRAY_HOST_DEVICE constexpr const std::array<scalar, 2>& mom_range()
            const {
            return m_mom_range;
        }
        DETRAY_HOST_DEVICE constexpr const point3& origin() const {
            return m_origin;
        }
        DETRAY_HOST_DEVICE constexpr const point3& origin_stddev() const {
            return m_origin_stddev;
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
                    ? point3{m_rnd_numbers.normal(m_cfg.origin()[0],
                                                  m_cfg.origin_stddev()[0]),
                             m_rnd_numbers.normal(m_cfg.origin()[1],
                                                  m_cfg.origin_stddev()[1]),
                             m_rnd_numbers.normal(m_cfg.origin()[2],
                                                  m_cfg.origin_stddev()[2])}
                    : m_cfg.origin();

            const std::array<scalar, 2>& phi_rng = m_cfg.phi_range();
            const std::array<scalar, 2>& theta_rng = m_cfg.theta_range();
            const std::array<scalar, 2>& mom_rng = m_cfg.mom_range();
            scalar p_mag{math::max(m_rnd_numbers({mom_rng[0], mom_rng[1]}),
                                   scalar{0.f})};
            scalar phi{std::clamp(m_rnd_numbers({phi_rng[0], phi_rng[1]}),
                                  phi_rng[0], phi_rng[1])};
            scalar theta{std::clamp(m_rnd_numbers({theta_rng[0], theta_rng[1]}),
                                    theta_rng[0], theta_rng[1])};
            const scalar sin_theta{math::sin(theta)};

            // Momentum direction from angles
            vector3 mom{math::cos(phi) * sin_theta, math::sin(phi) * sin_theta,
                        math::cos(theta)};
            // Magnitude of momentum
            vector::normalize(mom);

            mom = (m_cfg.is_pT() ? 1.f / sin_theta : 1.f) * p_mag * mom;

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
    /// If a concrete seed was given, use it otherwise prefer random seeds
    DETRAY_HOST_DEVICE
    constexpr random_track_generator(configuration cfg)
        : m_gen{(cfg.seed() != 0u) ? generator_t{cfg.seed()} : generator_t{}},
          m_cfg(cfg) {}

    /// Paramtetrized constructor for quick construction of simple tasks
    ///
    /// @note For more complex tasks, use the @c configuration type
    ///
    /// @param n_tracks the number of steps in the theta space
    /// @param mom_range the range of the track momentum (in GeV)
    /// @param charge charge of particle (e)
    DETRAY_HOST_DEVICE
    random_track_generator(
        std::size_t n_tracks,
        std::array<scalar, 2> mom_range = {1.f * unit<scalar>::GeV,
                                           1.f * unit<scalar>::GeV},
        scalar charge = -1.f * unit<scalar>::e)
        : m_gen{}, m_cfg{} {
        m_cfg.n_tracks(n_tracks);
        m_cfg.mom_range(mom_range);
        m_cfg.charge(charge);
    }

    /// Move constructor
    DETRAY_HOST_DEVICE
    random_track_generator(random_track_generator&& other)
        : m_gen(std::move(other.m_gen)), m_cfg(std::move(other.m_cfg)) {}

    /// Access the configuration
    DETRAY_HOST_DEVICE
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
