/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"

// System include(s).
#include <bits/stdc++.h>

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

}  // namespace detray