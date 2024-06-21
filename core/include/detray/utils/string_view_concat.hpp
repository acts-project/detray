/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <string>
#include <string_view>

namespace detray {
/// @brief Convenience class to statically concatenate two string views.
struct string_view_concat2 {
    std::string_view s1, s2;

    operator std::string() const { return std::string(s1) + std::string(s2); }
};
}  // namespace detray
