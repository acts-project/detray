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

/// @brief A proto landmark class as a simple translation layer from a
/// description of a point.
template <concepts::point3D point3_t>
struct trajectory {
    std::vector<point3_t> _points;
    std::string _name;
    actsvg::style::stroke _stroke = actsvg::style::stroke();
};
}  // namespace detray::svgtools::meta::proto
