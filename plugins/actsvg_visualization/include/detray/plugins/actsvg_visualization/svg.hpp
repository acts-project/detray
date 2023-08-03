#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/surface.hpp"
#include "detray/plugins/actsvg_visualization/volume.hpp"
#include "detray/plugins/actsvg_visualization/detector.hpp"
#include "detray/plugins/actsvg_visualization/options.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/detector_writer.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization {

/// @brief Convert a detector or its components to SVG.
template <typename detector_t>
class detector_visualizer {

    public:

    detector_visualizer() = delete;

    constexpr detector_visualizer(const detector_t& det)
        : _detector{det} {}

    /// @brief Converts a detray surface in the detector to an actsvg svg.
    /// @param context the geometry context.
    /// @param options the visualization options.
    /// @note Does not override visualization options.
    /// @returns SVG of the detector's surface in x-y view.
    inline auto xy_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& surface, const detector::detector_options& options = default_options) const{
        return surface::to_svg(context, _xy_view, surface, options.v_options.s_options, options.v_options.p_options, "");
    }

    /// @brief Converts a detray volume in the detector to an actsvg svg.
    /// @param context the geometry context.
    /// @param options the visualization options.
    /// @note Does not override visualization options.
    /// @returns SVG of the detector's volume in x-y view.
    inline auto xy_volume(const typename detector_t::geometry_context& context, const detray::detector_volume<detector_t>& volume, const detector::detector_options& options = default_options) const{
        return volume::to_svg(context, _xy_view, _detector, volume, options.v_options, "");
    }

    /// @brief Converts a detray detector to an actsvg svg.
    /// @param context the geometry context.
    /// @param options the visualization options.
    /// @returns SVG of the detector in x-y view.
    inline auto xy(const typename detector_t::geometry_context& context, const detector::detector_options& options = default_options) const{
        return detector::to_svg(context, _xy_view, _detector, options, "");
    }

    /// @brief Converts a detray surface in the detector to an actsvg svg.
    /// @param context the geometry context.
    /// @param options the visualization options.
    /// @note Does not override visualization options.
    /// @returns SVG of the detector's surface in z-r view.
    inline auto zr_surface(const typename detector_t::geometry_context& context, const detray::surface<detector_t>& surface, const detector::detector_options& options = default_options) const{
        return surface::to_svg(context, _zr_view, surface, options.v_options.s_options, options.v_options.p_options, "");
    }

    /// @brief Converts a detray volume in the detector to an actsvg svg.
    /// @param context the geometry context.
    /// @param options the visualization options.
    /// @note Does not override visualization options.
    /// @returns SVG of the detector's volume in z-r view.
    inline auto zr_volume(const typename detector_t::geometry_context& context, const detray::detector_volume<detector_t>& volume, const detector::detector_options& options = default_options) const{
        return volume::to_svg(context, _zr_view, _detector, volume, options.v_options, "");
    }

    /// @brief Converts a detray detector to an actsvg svg.
    /// @param context the geometry context.
    /// @param options the visualization options.
    /// @return SVG of the detector in z-r view.
    inline auto zr(const typename detector_t::geometry_context& context, const detector::detector_options& options = default_options) const{
        return detector::to_svg(context, _zr_view, _detector, options, "");
    }

    private:

    // x-y view.
    const actsvg::views::x_y _xy_view;
    // z-r view.
    const actsvg::views::z_r _zr_view;
    /// Access to the detector stores.
    const detector_t& _detector;

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
inline void write_svg(const std::string& path, const std::initializer_list<actsvg::svg::object>& svgs){
    std::vector<actsvg::svg::object> arg = svgs;
    write_svg(path, arg);
}

/// @brief Writes an svg objects to a file.
inline void write_svg(const std::string& path, const actsvg::svg::object& svg){
    write_svg(path, std::array{svg});
}

}  // namespace detray::actsvg_visualization