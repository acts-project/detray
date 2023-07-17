#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg/portal_conversion.hpp"
#include "detray/plugins/actsvg/surface_conversion.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization {

    const visualization_style default_style;

    struct visualization_style
    {
        actsvg::style::color portal_color = {{178, 102, 255}, 0.5};
        actsvg::style::color surface_color = {{102, 102, 255}, 0.5};
    }

    template <typename detector_t, typename view_t>
    actsvg::svg::object svg(
    const std::string object_name
    const detray::surface<detector_t>& d_surface,
    const typename detector_t::geometry_context& context,
    const view_t view,
    const visualization_style style = default_style
    )
    const 
    {
        if (d_surface.is_portal()){
            return actsvg::display::portal(object_name, convert_portal(d_surface, context), view);
        }
        return actsvg::display::surface(object_name, convert_surface(d_surface, context), view);
    }
}  // namespace detray::actsvg_visualization