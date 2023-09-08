/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <random>
#include <type_traits>

namespace detray::detail {

/// Wrapper for random number generatrion
template <typename scalar_t = detray::scalar,
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

}  // namespace detray::detail
