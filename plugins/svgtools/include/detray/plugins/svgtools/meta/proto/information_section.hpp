/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <vector>
#include <string>

namespace detray::svgtools::meta::proto {

/// @brief A proto intersection record class as a simple translation layer from a intersection record description.
template <typename point3_t>
class information_section{
    public:
    std::vector<std::string> _info;
    std::string _title;
    std::array<int, 3> _color;
    point3_t _position;
};
}