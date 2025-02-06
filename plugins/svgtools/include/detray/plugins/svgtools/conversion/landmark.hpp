// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Project include(s)
#include "detray/plugins/svgtools/meta/proto/landmark.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"

namespace detray::svgtools::conversion {

/// @returns The proto landmark of a detray point.
template <concepts::point3D point3_t>
inline auto landmark(const point3_t& position,
                     const styling::landmark_style& style =
                         styling::svg_default::landmark_style) {

    using p_landmark_t = svgtools::meta::proto::landmark<point3_t>;

    p_landmark_t p_lm;
    p_lm._position = position;

    svgtools::styling::apply_style(p_lm, style);

    return p_lm;
}

}  // namespace detray::svgtools::conversion
