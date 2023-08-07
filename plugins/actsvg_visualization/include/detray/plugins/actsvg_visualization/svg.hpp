#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/volume.hpp"
#include "detray/plugins/actsvg_visualization/styling/styling.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/detector_writer.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization {

/// @brief Conversion of detector or volumes/surfaces to SVG.
template <typename detector_t>
class svg_converter{

    public:

    svg_converter() = delete;

    svg_converter(const detector_t& detector, const typename detector_t::name_map& name_map) 
    : _detector{detector}, _name_map{name_map} {}

    svg_converter(const detector_t& detector, const typename detector_t::name_map& name_map, const styling::style& style) 
    : _detector{detector}, _name_map{name_map}, _style{style} {}

    /// @brief Converts a detray surface in the detector to an actsvg svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param index the index of the surface in the detector.
    /// @returns SVG of the detector's surface in a x-y view.
    inline auto xy_surface(const std::string& identification, const typename detector_t::geometry_context& context, const size_t index) const {
        return svg_surface(identification, context, index, _xy_view);
    }

    /// @brief Converts a collection of detray surfaces in the detector to an actsvg svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param indices the collection of surface indices in the detector to convert.
    /// @returns SVG of the given detector surfaces in a x-y view.
    template <typename iterator_t>
    inline auto xy_surfaces(const std::string& identification, const typename detector_t::geometry_context& context, const iterator_t& indices) const {
        return svg_surfaces(identification, context, indices, _xy_view);
    }

    /// @brief Converts a detray volume in the detector to an actsvg svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param index the index of the volume in the detector.
    /// @returns SVG of the detector's volume in a x-y view.
    inline auto xy_volume(const std::string& identification, const typename detector_t::geometry_context& context, const size_t index) const {
        return svg_volume(identification, context, index, _xy_view);
    }

    /// @brief Converts a collection of detray volumes in the detector to an actsvg svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param indices the collection of volume indices in the detector to convert.
    /// @returns SVG of the given detector volumes in a x-y view.
    template <typename iterator_t>
    inline auto xy_volumes(const std::string& identification, const typename detector_t::geometry_context& context, const iterator_t& indices) const {
        return svg_volumes(identification, context, indices, _xy_view);
    }

    /// @brief Converts a detray surface detector to an actsvg svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @returns SVG of the detector in a x-y view.
    inline auto xy_detector(const std::string& identification, const typename detector_t::geometry_context& context) const {
        return svg(identification, context, _xy_view);
    }

    /// @brief Converts a detray surface in the detector to an actsvg svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param index the index of the surface in the detector.
    /// @returns SVG of the detector's surface in a z-r view.
    inline auto zr_surface(const std::string& identification, const typename detector_t::geometry_context& context, const size_t index) const {
        return svg_surface(identification, context, index, _zr_view);
    }

    /// @brief Converts a collection of detray surfaces in the detector to an actsvg svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param indices the collection of surface indices in the detector to convert.
    /// @returns SVG of the given detector surfaces in a z-r view.
    template <typename iterator_t>
    inline auto zr_surfaces(const std::string& identification, const typename detector_t::geometry_context& context, const iterator_t& indices) const {
        return svg_surfaces(identification, context, indices, _zr_view);
    }

    /// @brief Converts a detray volume in the detector to an actsvg svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param index the index of the volume in the detector.
    /// @returns SVG of the detector's volume in a z-r view.
    inline auto zr_volume(const std::string& identification, const typename detector_t::geometry_context& context, const size_t index) const {
        return svg_volume(identification, context, index, _zr_view);
    }

    /// @brief Converts a collection of detray volumes in the detector to an actsvg svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param indices the collection of volume indices in the detector to convert.
    /// @returns SVG of the given detector volumes in a z-r view.
    template <typename iterator_t>
    inline auto zr_volumes(const std::string& identification, const typename detector_t::geometry_context& context, const iterator_t& indices) const {
        return svg_volumes(identification, context, indices, _zr_view);
    }

    /// @brief Converts a detray surface detector to an actsvg svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @returns SVG of the detector in a z-r view.
    inline auto zr_detector(const std::string& identification, const typename detector_t::geometry_context& context) const {
        return svg(identification, context, _zr_view);
    }

    private:

    template <typename view_t>
    inline auto svg_surface(const std::string& identification, const typename detector_t::geometry_context& context, const size_t index, const view_t& view) const {
        const auto surface = _detector.surface_by_index(static_cast<detray::dindex>(index));
        if (surface.is_portal())
        {
            auto p_portal = proto::portal(context, _detector, surface);
            styling::apply_style(p_portal, _style);
            return actsvg::display::portal(identification, p_portal, view);
        }
        auto p_surface = proto::surface(context, surface);
        styling::apply_style(p_surface, _style);
        return actsvg::display::surface(identification, p_surface, view);
    }

    template <typename iterator_t, typename view_t>
    inline auto svg_surfaces(const std::string& identification, const typename detector_t::geometry_context& context, const iterator_t& indices, const view_t& view) const {
        actsvg::svg::object ret;
        ret._tag = "g";
        ret._id = identification;
        for (const auto index : indices){
            const auto svg = svg_surface(identification + "_surface" + std::to_string(index), context, index, view);
            ret.add_object(svg);
        }
        return ret;
    }

    template <typename view_t>
    inline auto svg_volume(const std::string& identification, const typename detector_t::geometry_context& context, const size_t index, const view_t& view) const {
        const auto volume = _detector.volume_by_index(static_cast<detray::dindex>(index));
        auto p_volume = proto::volume(context, _detector, volume);
        styling::apply_style(p_volume, _style);
        return actsvg::display::volume(identification, p_volume, view);
    }

    template <typename iterator_t, typename view_t>
    inline auto svg_volumes(const std::string& identification, const typename detector_t::geometry_context& context, const iterator_t& indices, const view_t& view) const {
        actsvg::svg::object ret;
        ret._tag = "g";
        ret._id = identification;
        for (const auto index : indices){
            const auto svg = svg_volume(identification + "_volume" + std::to_string(index), context, index, view);
            ret.add_object(svg);
        }
        return ret;
    }

    template <typename view_t>
    inline auto svg(const std::string& identification, const typename detector_t::geometry_context& context, const view_t& view) const {
        actsvg::svg::object ret;
        ret._tag = "g";
        ret._id = identification;
        for (size_t index = 0; index < _detector.volumes().size(); index++){
            const auto svg = svg_volume(identification + "_volume" + std::to_string(index), context, index, view);
            ret.add_object(svg);
        }
        return ret;
    }

    const detector_t& _detector;
    const typename detector_t::name_map& _name_map;
    const styling::style _style = styling::default_style;
    const actsvg::views::x_y _xy_view;
    const actsvg::views::z_r _zr_view;
};

/// @brief Writes a collection of svgs objects to a single file.
template <typename container_t>
inline void write_svg(const std::string& path, const container_t& svgs){
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
/// @note To avoid conflict, the ids of the svg objects must be unique.
inline void write_svg(const std::string& path, const std::initializer_list<actsvg::svg::object>& svgs){
    std::vector<actsvg::svg::object> arg = svgs;
    write_svg(path, arg);
}

/// @brief Writes an svg object to a file.
/// @note To avoid conflict, the ids of the svg objects must be unique.
inline void write_svg(const std::string& path, const actsvg::svg::object& svg){
    write_svg(path, std::array{svg});
}

}  // namespace detray::actsvg_visualization