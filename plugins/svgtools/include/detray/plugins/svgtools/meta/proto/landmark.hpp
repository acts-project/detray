#pragma once

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <string>

namespace detray::svgtools::meta::proto {

template <typename point3_t>
class landmark{
    public:
    point3_t _position;
    std::string _name;
    actsvg::style::marker _marker{"x", 1.};
};
}  