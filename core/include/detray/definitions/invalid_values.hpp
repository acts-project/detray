/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <limits>

#include "detray/definitions/detray_qualifiers.hpp"

template <typename T>
DETRAY_HOST_DEVICE inline T invalid_value() {
    return T::invalid_value();
}

// specialization for int
template <>
DETRAY_HOST_DEVICE inline int invalid_value() {
    return std::numeric_limits<int>::max();
}

// specialization for unsigned int
template <>
DETRAY_HOST_DEVICE inline unsigned int invalid_value() {
    return std::numeric_limits<unsigned int>::max();
}

// specialization for long unsigned int
template <>
DETRAY_HOST_DEVICE inline long unsigned int invalid_value() {
    return std::numeric_limits<long unsigned int>::max();
}
