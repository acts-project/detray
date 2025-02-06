// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <string>
#include <vector>

namespace detray::svgtools::meta::proto {

/// @brief A proto intersection record class as a simple translation layer from
/// a intersection record description.
template <concepts::point3D point3_t>
struct information_section {
    std::vector<std::string> _info;
    std::string _title;
    std::array<int, 3> _color;
    point3_t _position;
};

}  // namespace detray::svgtools::meta::proto
