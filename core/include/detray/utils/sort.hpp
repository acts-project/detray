/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <algorithm>
#include <functional>

namespace detray {

template <class RandomIt>
DETRAY_HOST_DEVICE inline void insertion_sort(RandomIt first, RandomIt last) {
    for (auto it = first; it != last; it++) {
        // Searching the upper bound, i.e., first
        // element greater than *it from beginning
        auto const insertion_point = std::upper_bound(first, it, *it);

        // Shifting the unsorted part
        std::rotate(insertion_point, it, it + 1);
    }
}

// Function to sort the array
template <template <typename...> class vector_t, typename TYPE>
DETRAY_HOST_DEVICE inline void insertion_sort(vector_t<TYPE> &vec) {
    insertion_sort(vec.begin(), vec.end());
}

template <class RandomIt, class Comp = std::less<void>>
DETRAY_HOST_DEVICE inline void selection_sort(RandomIt first, RandomIt last,
                                              Comp &&comp = Comp()) {
    for (RandomIt i = first; i < (last - 1); ++i) {
        RandomIt k = i;

        for (RandomIt j = i + 1; j < last; ++j) {
            if (comp(*j, *k)) {
                k = j;
            }
        }

        if (k != i) {
            auto t = *i;
            *i = *k;
            *k = t;
        }
    }
}

// Function to sort the array
template <template <typename...> class vector_t, typename TYPE>
DETRAY_HOST_DEVICE inline void selection_sort(vector_t<TYPE> &vec) {
    selection_sort(vec.begin(), vec.end());
}

}  // namespace detray
