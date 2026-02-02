/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <cstdint>
#include <ostream>

namespace detray::actor {

enum class status : std::uint8_t {
    e_notify = 0u,
    e_success = 1u,
    e_unknown = 2u,
    e_failure = 3u,
};

// Print the values of an enum by identifier
#define ENUM_PRINT(x) \
    case x:           \
        os << #x;     \
        break

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, status s) {

    switch (s) {
        using enum status;
        ENUM_PRINT(e_notify);
        ENUM_PRINT(e_success);
        ENUM_PRINT(e_unknown);
        ENUM_PRINT(e_failure);
    }
    return os;
}

#undef ENUM_PRINT
}  // namespace detray::actor
