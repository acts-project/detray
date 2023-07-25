#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg/portal_conversion.hpp"
#include "detray/plugins/actsvg/surface_conversion.hpp"
#include "detray/plugins/actsvg/volume_conversion.hpp"
#include "detray/plugins/actsvg/detector_conversion.hpp"
#include "detray/plugins/actsvg/proto_styling.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/detector_writer.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization {

template <typename shape_t, typename view_t>
actsvg::svg::object svg(
const std::string& object_name,
const detray::mask<shape_t>& mask,
const view_t& view,
const detector_style& style = default_detector_style
)
{
    auto p_surface = convert_mask(mask);
    apply_style(p_surface, style);
    return actsvg::display::surface(object_name, p_surface, view);
}

/// @brief Converts the detray surface/portal to svg.
/// @returns svg.
template <typename detector_t, typename view_t>
actsvg::svg::object svg(
const std::string& object_name,
const detector_t detector,
const detray::surface<detector_t>& d_surface,
const typename detector_t::geometry_context& context,
const view_t& view,
const detector_style& style = default_detector_style
)
{
    if (d_surface.is_portal()){
        auto p_portal = convert_portal(detector, d_surface, context);
        apply_style(p_portal, style);
        return actsvg::display::portal(object_name, p_portal, view);
    }

    auto p_surface = convert_surface(d_surface, context);
    apply_style(p_surface, style);
    return actsvg::display::surface(object_name, p_surface, view);
}

/// @brief Converts the detray volume to svg.
/// @returns svg.
template <typename detector_t, typename view_t>
actsvg::svg::object svg(
const std::string& object_name,
const detector_t detector,
const detray::detector_volume<detector_t>& d_volume,
const typename detector_t::geometry_context& context,
const view_t& view,
const detector_style& style = default_detector_style
)
{
    auto p_volume = convert_volume(detector, d_volume, context);
    apply_style(p_volume, style);
    return actsvg::display::volume(object_name, p_volume, view);
}

/// @brief Converts the detray detector to svg.
/// @returns svg.
template <typename detector_t, typename view_t>
actsvg::svg::object svg(
const std::string object_name,
const detector_t detector,
const typename detector_t::geometry_context& context,
const view_t view,
const detector_style& style = default_detector_style
)
{
    auto p_detector = convert_detector(detector, context);
    apply_style(p_detector, style);
    return actsvg::display::detector(object_name, p_detector, view);
}

/// @brief Writes a collection of svgs objects to a single file.
template <typename container_t>
void write_svg(const container_t& svgs, const std::string& path){
    actsvg::svg::file file;
    for (const actsvg::svg::object& obj : svgs){
        file.add_object(obj);
    }
    detray::io::detail::file_handle stream{path,
                                        "",
                                        std::ios::out | std::ios::trunc};
    *stream << file;
}

/// @brief Writes an svg objects to a file.
void write_svg(const actsvg::svg::object& svg, const std::string& path){
    write_svg(std::array{svg}, path);
}

}  // namespace detray::actsvg_visualization