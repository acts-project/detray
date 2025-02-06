// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Project inlude(s)
#include "detray/plugins/svgtools/meta/proto/eta_lines.hpp"

// Actsvg include(s)
#include "actsvg/display/geometry.hpp"

// System include(s)
#include <tuple>

namespace detray::svgtools::meta::display {

/// @brief Converts a proto eta_line to an SVG object.
inline auto eta_lines(const std::string& id,
                      const svgtools::meta::proto::eta_lines& el) {

    return actsvg::display::eta_lines(
        id, el._z, el._r,
        {std::make_tuple(el._values_main, el._stroke_main, el._show_label,
                         el._label_font),
         std::make_tuple(el._values_half, el._stroke_half, false,
                         el._label_font)});
}

}  // namespace detray::svgtools::meta::display
