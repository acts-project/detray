// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <string>

namespace detray::svgtools::meta::proto {

/// @brief A proto landmark class as a simple translation layer from a
/// description of a point.
template <concepts::point3D point3_t>
struct landmark {
    point3_t _position{0.f, 0.f, 0.f};
    std::string _name{"unknown landmark"};
    actsvg::style::marker _marker{"x", 1.f};
};

}  // namespace detray::svgtools::meta::proto
