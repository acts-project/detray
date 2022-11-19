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
        : stddev(smearer.stddev), generator(smearer.generator) {}

    void set_seed(const std::size_t sd) { generator.seed(sd); }

    std::array<scalar_t, 2> stddev;
    std::random_device rd{};
    std::mt19937 generator{rd()};

    std::array<scalar_t, 2> get_offset() {
        return {std::normal_distribution<scalar_t>(0, stddev[0])(generator),
                std::normal_distribution<scalar_t>(0, stddev[1])(generator)};
    }

    template <typename local_parameter_t>
    std::array<scalar_t, 2> operator()(const std::string& name,
                                       const std::size_t& meas_dim,
                                       const std::array<scalar_t, 2>& offset,
                                       const local_parameter_t& param) {

        std::array<scalar_t, 2> ret;

        const auto local = param.local();

        // The radial element of line measurement (ret[0]) always have a
        // positive value
        if (name == "line") {
            ret[0] = std::abs(local[0]) + offset[0];
            if (ret[0] < 0.) {
                ret[0] = 0.;
            }
        } else {
            ret[0] = local[0] + offset[0];
        }

        // If the measurement dimension is 1, the second element is null
        if (meas_dim == 1) {
            ret[1] = 0.;
        } else if (meas_dim == 2) {
            ret[1] = local[1] + offset[1];
        }

        return ret;
    }
};

}  // namespace detray