/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <algorithm>
#include <numeric>

namespace detray {

template <typename scalar_t, template <typename...> class vector_t>
scalar_t get_variance(const vector_t<scalar_t>& v) {
    scalar_t sum = std::accumulate(v.begin(), v.end(), 0.0);
    scalar_t mean = sum / v.size();
    std::vector<scalar_t> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),
                   [mean](scalar_t x) { return x - mean; });
    scalar_t sq_sum =
        std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    scalar_t variance = sq_sum / v.size();
    return variance;
}

}  // namespace detray