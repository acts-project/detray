/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <ratio>
#include <utility>

namespace detray {

template <typename... ratios>
struct ratio_sum_helper;

template <typename ratio>
struct ratio_sum_helper<ratio> {
    using first_two_sum = ratio;
};

template <typename ratio1, typename ratio2, typename... ratios>
struct ratio_sum_helper<ratio1, ratio2, ratios...> {

    using first_two_sum = std::ratio_add<ratio1, ratio2>;

    using next_helper = ratio_sum_helper<first_two_sum, ratios...>;

    static constexpr bool is_done = (sizeof...(ratios) == 0);

    using sum =
        typename std::conditional_t<is_done, first_two_sum,
                                    typename next_helper::first_two_sum>;
};

template <typename... ratios>
struct ratio_sum {

    using helper = ratio_sum_helper<ratios...>;

    using sum = typename helper::sum;

    static constexpr bool is_sum_one = (sum::num == sum::den);
};

}  // namespace detray