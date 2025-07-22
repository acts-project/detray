/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <array>
#include <limits>
#include <random>

namespace detray::detail {

/// Wrapper for CPU random number generatrion for the @c random_track_generator
template <typename scalar_t = scalar,
          typename distribution_t = std::uniform_real_distribution<scalar_t>,
          typename engine_t = std::mt19937_64>
struct random_numbers {

    using distribution_type = distribution_t;
    using engine_type = engine_t;
    using seed_type = typename engine_t::result_type;

    std::seed_seq m_seeds;
    engine_t m_engine;

    /// Default seed
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST
    random_numbers()
        : m_seeds{random_numbers::default_seed()}, m_engine{m_seeds} {}
#endif // DETRAY_COMPILE_VITIS

    /// Different seed @param s for every instance
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST
    random_numbers(seed_type s) : m_seeds{s}, m_engine{m_seeds} {}
#endif // DETRAY_COMPILE_VITIS

    /// More entropy in seeds from collection @param s
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST
    random_numbers(const std::vector<seed_type>& s)
        : m_seeds{s.begin(), s.end()}, m_engine{m_seeds} {}
#endif // DETRAY_COMPILE_VITIS

    /// Copy constructor
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST
    random_numbers(random_numbers&& other)
        : m_engine(std::move(other.m_engine)) {}
#endif // DETRAY_COMPILE_VITIS

    /// Generate random numbers in a given range
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST auto operator()(const std::array<scalar_t, 2> range = {
                                    -std::numeric_limits<scalar_t>::max(),
                                    std::numeric_limits<scalar_t>::max()}) {
        const scalar_t min{range[0]}, max{range[1]};
        assert(min <= max);

        // Uniform
        if constexpr (std::is_same<
                          distribution_t,
                          std::uniform_real_distribution<scalar_t>>::value ) {
            return distribution_t(min, max)(m_engine);

            // Normal
        } else if constexpr (std::is_same<
                                 distribution_t,
                                 std::normal_distribution<scalar_t>>::value ) {
            scalar_t mu{min + 0.5f * (max - min)};
            return distribution_t(mu, 0.5f / 3.0f * (max - min))(m_engine);
        }
    }
#endif // DETRAY_COMPILE_VITIS

    /// Explicit normal distribution around a @param mean and @param stddev
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST auto normal(const scalar_t mean, const scalar_t stddev) {
        return std::normal_distribution<scalar_t>(mean, stddev)(m_engine);
    }
#endif // DETRAY_COMPILE_VITIS

    /// 50:50 coin toss
#ifndef DETRAY_COMPILE_VITIS
    DETRAY_HOST std::uint8_t coin_toss() {
        return std::uniform_int_distribution<std::uint8_t>(0u, 1u)(m_engine);
    }
#endif // DETRAY_COMPILE_VITIS

    /// Get the default seed of the engine
    static constexpr seed_type default_seed() { return engine_t::default_seed; }
};

}  // namespace detray::detail
