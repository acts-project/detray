/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/tracks/bound_track_parameters.hpp"

// System include(s).
#include <array>
#include <random>
#include <string>

namespace detray {

template <typename transform3_t>
struct measurement_smearer {

    using matrix_operator = typename transform3_t::matrix_actor;
    using scalar_type = typename transform3_t::scalar_type;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;

    measurement_smearer(const scalar_type stddev_local0,
                        const scalar_type stddev_local1)
        : stddev({stddev_local0, stddev_local1}) {}

    measurement_smearer(measurement_smearer& smearer)
        : stddev(smearer.stddev), generator(smearer.generator) {}

    void set_seed(const uint_fast64_t sd) { generator.seed(sd); }

    std::array<scalar_type, 2> stddev;
    std::random_device rd{};
    std::mt19937_64 generator{rd()};

    std::array<scalar_type, 2> get_offset() {
        return {
            std::normal_distribution<scalar_type>(0.f, stddev[0])(generator),
            std::normal_distribution<scalar_type>(0.f, stddev[1])(generator)};
    }

    template <typename mask_t>
    std::array<scalar_type, 2> operator()(
        const mask_t& mask, const std::array<scalar_type, 2>& offset,
        const bound_track_parameters<transform3_t>& bound_params) {

        std::array<scalar_type, 2> ret{0.f, 0.f};

        matrix_type<mask_t::shape::meas_dim, 1u> meas =
            mask.projection_matrix(bound_params) * bound_params.vector();

        if constexpr (mask_t::shape::meas_dim == 1u) {
            ret[0u] = matrix_operator().element(meas, 0u, 0u) + offset[0];
            ret[1u] = scalar_type{0.f};
            // Special treatment for line coordinate: we don't allow negative
            // radial distance
            if (mask_t::shape::name == "line" &&
                mask_t::shape::normal_order == true) {
                ret[0u] = std::max(ret[0u], static_cast<scalar_type>(0.f));
            }

        } else if (mask_t::shape::meas_dim == 2u) {
            ret[0u] = matrix_operator().element(meas, 0u, 0u) + offset[0];
            ret[1u] = matrix_operator().element(meas, 1u, 0u) + offset[1];

            // Special treatment for line coordinate: we don't allow negative
            // radial distance
            if (mask_t::shape::name == "line") {

                if constexpr (mask_t::shape::normal_order == true) {
                    ret[0u] = std::max(ret[0u], static_cast<scalar_type>(0.f));
                } else {
                    ret[1u] = std::max(ret[1u], static_cast<scalar_type>(0.f));
                }
            }
        }

        return ret;
    }
};

}  // namespace detray