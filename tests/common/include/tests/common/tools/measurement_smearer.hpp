/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <random>

namespace detray {

template <typename scalar_t>
struct measurement_smearer {

    measurement_smearer(scalar_t stddev_local0, scalar_t stddev_local1)
        : stddev({stddev_local0, stddev_local1}) {}

    measurement_smearer(measurement_smearer<scalar_t>& smearer)
        : stddev(smearer.stddev) {}

    std::array<scalar_t, 2> stddev;

    std::random_device rd{};
    std::mt19937 generator{rd()};

    template <int ID>
    scalar_t get() {
        static_assert(ID == 0 || ID == 1);

        if constexpr (ID == 0) {
            return std::normal_distribution<scalar_t>(0, stddev[0])(generator);
        } else if (ID == 1) {
            return std::normal_distribution<scalar_t>(0, stddev[1])(generator);
        }
    }
};

}  // namespace detray