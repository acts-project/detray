// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <string>
#include <string_view>

namespace detray {
/// @brief Convenience class to statically concatenate two string views.
struct string_view_concat2 {
    std::string_view s1;
    std::string_view s2;

    explicit operator std::string() const {
        return std::string(s1) + std::string(s2);
    }
};
}  // namespace detray
