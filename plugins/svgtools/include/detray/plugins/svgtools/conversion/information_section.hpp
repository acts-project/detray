/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once 

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/svgtools/meta/proto/information_section.hpp"

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "actsvg/meta.hpp"

// System include(s)
#include <array>
#include <string>
#include <iostream>

namespace detray::svgtools::conversion {

template <typename point3_t, typename detector_t>
inline auto information_section(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& d_surface)
{
    const auto title = "My Title";
    const std::vector<std::string> info = {"Hello", "World"};
    svgtools::meta::proto::information_section<point3_t> is;
    is._title = title;
    is._info = info;
    is._position = {0, 0, 0};
    return is;
}


}