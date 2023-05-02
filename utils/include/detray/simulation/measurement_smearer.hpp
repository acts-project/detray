/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <array>
#include <random>
#include <string>

namespace detray {

template <typename scalar_t>
struct measurement_smearer {

    measurement_smearer(scalar_t stddev_local0, scalar_t stddev_local1)
        : stddev({stddev_local0, stddev_local1}) {}

    measurement_smearer(measurement_smearer<scalar_t>& smearer)
        : stddev(smearer.stddev), generator(smearer.generator) {}

    void set_seed(const uint_fast64_t sd) { generator.seed(sd); }

    std::array<scalar_t, 2> stddev;
    std::random_device rd{};
    std::mt19937_64 generator{rd()};

    std::array<scalar_t, 2> get_offset() {
        return {std::normal_distribution<scalar_t>(0.f, stddev[0])(generator),
                std::normal_distribution<scalar_t>(0.f, stddev[1])(generator)};
    }

    template <typename local_parameter_t>
    std::array<scalar_t, 2> operator()(const std::string& name,
                                       const std::size_t& meas_dim,
                                       const std::array<scalar_t, 2>& offset,
                                       const local_parameter_t& param) {

        std::array<scalar_t, 2> ret;

        const auto local = param.bound_local();

        // The radial element of line measurement (ret[0]) always have a
        // positive value
        if (name == "line") {
            ret[0] = std::abs(local[0]) + offset[0];
            if (ret[0] < 0.f) {
                ret[0] = 0.f;
            }
        } else {
            ret[0] = local[0] + offset[0];
        }

        // If the measurement dimension is 1, the second element is null
        if (meas_dim == 1u) {
            ret[1] = 0.f;
        } else if (meas_dim == 2u) {
            ret[1] = local[1] + offset[1];
        }

        return ret;
    }
};

}  // namespace detray