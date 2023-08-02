#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/portal.hpp"
#include "detray/plugins/actsvg_visualization/surface.hpp"
#include "detray/plugins/actsvg_visualization/volume.hpp"
#include "detray/plugins/actsvg_visualization/detector.hpp"
#include "detray/plugins/actsvg_visualization/proto_styling.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/detector_writer.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <assert.h>
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization {

inline void assert_legal_name(const std::string& name){
    auto legal_name = name;
    legal_name.erase(remove_if(legal_name.begin(), legal_name.end(), isspace), legal_name.end());
    assert(name == legal_name);
}

/// @brief Converts the detray surface/portal to svg.
/// @returns svg.
template <typename detector_t, typename view_t>
actsvg::svg::object to_svg(
const typename detector_t::geometry_context& context,
const detector_t detector,
const detray::surface<detector_t>& d_surface,
const view_t& view,
const detector_style& style = default_detector_style
)
{
    const auto object_name = "surface";
    if (d_surface.is_portal()){
        auto p_portal = portal::to_proto(context, d_surface);
        apply_style(p_portal, style);
        return actsvg::display::portal(object_name, p_portal, view);
    }

    auto p_surface = surface::to_proto(context, d_surface);
    apply_style(p_surface, style);
    return actsvg::display::surface(object_name, p_surface, view);
}

/// @brief Converts the detray volume to svg.
/// @returns svg.
template <typename detector_t, typename view_t>
actsvg::svg::object to_svg(
const typename detector_t::geometry_context& context,
const detector_t detector,
const detray::detector_volume<detector_t>& d_volume,
const view_t& view,
const detector_style& style = default_detector_style
)
{
    const auto object_name = "volume";
    auto p_volume = volume::to_proto(context, detector, d_volume);
    apply_style(p_volume, style);
    return actsvg::display::volume(object_name, p_volume, view);
}

/// @brief Converts the detray detector to svg.
/// @returns svg.
template <typename detector_t, typename view_t>
actsvg::svg::object to_svg(
const typename detector_t::geometry_context& context,
const detector_t detector,
const view_t view,
const detector_style& style = default_detector_style
)
{
    const auto object_name = "detector";
    auto p_detector = detector::to_proto(context, detector);
    apply_style(p_detector, style);
    return actsvg::display::detector(object_name, p_detector, view);
}

/// @brief Writes a collection of svgs objects to a single file.
template <typename container_t>
void write_svg(const std::string& path, const container_t& svgs){
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
void write_svg(const std::string& path, const std::initializer_list<actsvg::svg::object>& svgs){
    std::vector<actsvg::svg::object> arg = svgs;
    write_svg(path, arg);
}

/// @brief Writes an svg objects to a file.
void write_svg(const std::string& path, const actsvg::svg::object& svg){
    write_svg(path, std::array{svg});
}

}  // namespace detray::actsvg_visualization