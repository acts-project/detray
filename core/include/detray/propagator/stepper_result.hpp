/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/tracks/free_track_parameters.hpp"
#include "detray/utils/invalid_values.hpp"

namespace detray::stepping {

template <concepts::algebra algebra_t, std::size_t N = 1u>
struct alignas(16) result {
    using scalar_t = dscalar<algebra_t>;

    /// TODO: replace with @c dsimd<algebra_t, T>
    template <typename T>
    using simd_t = std::array<T, N>;

    DETRAY_HOST_DEVICE
    free_track_parameters<algebra_t> free_param(
        const std::size_t i = 0u) const {
        return free_track_parameters<algebra_t>{
            free_vector<algebra_t>{pos0[i], pos1[i], pos2[i], time[i], dir0[i],
                                   dir1[i], dir2[i], qop[i]}};
    };

    DETRAY_HOST_DEVICE
    void free_param(const free_track_parameters<algebra_t>& param,
                    const std::size_t i = 0u) {
        assert(i < N);

        pos0[i] = param[0u];
        pos1[i] = param[1u];
        pos2[i] = param[2u];
        time[i] = param[3u];
        dir0[i] = param[4u];
        dir1[i] = param[5u];
        dir2[i] = param[6u];
        qop[i] = param[7u];
    };

    // private:
    simd_t<scalar_t> pos0{};
    simd_t<scalar_t> pos1{};
    simd_t<scalar_t> pos2{};

    simd_t<scalar_t> time{};

    simd_t<scalar_t> dir0{};
    simd_t<scalar_t> dir1{};
    simd_t<scalar_t> dir2{};

    simd_t<scalar_t> qop{};
};

}  // namespace detray::stepping
