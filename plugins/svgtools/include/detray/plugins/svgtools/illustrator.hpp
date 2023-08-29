/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/svgtools/conversion/grid.hpp"
#include "detray/plugins/svgtools/conversion/information_section.hpp"
#include "detray/plugins/svgtools/conversion/intersection_record.hpp"
#include "detray/plugins/svgtools/conversion/landmark.hpp"
#include "detray/plugins/svgtools/conversion/surface.hpp"
#include "detray/plugins/svgtools/conversion/trajectory.hpp"
#include "detray/plugins/svgtools/conversion/volume.hpp"
#include "detray/plugins/svgtools/meta/display/geometry.hpp"
#include "detray/plugins/svgtools/meta/display/information.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"
#include "detray/plugins/svgtools/utils/volume_utils.hpp"
#include "detray/utils/ranges.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"

// System include(s)
#include <array>
#include <vector>

namespace detray::svgtools {

/// @brief SVG generator for a detector and related entities.
/// @note To view information boxes, they must be enabled in the constructor.
/// Furthermore the svg viewer (opening the file after it is created) must support animations.
template <typename detector_t>
class illustrator {

    public:
    illustrator() = delete;

    /// @param detector the detector
    /// @param context the geometry context
    /// @note information boxes are disabled unless expicitly enabled (use another constructor)
    illustrator(const detector_t& detector,
                const geometry_context& context)
        : _detector{detector}, _context{context} {}

    /// @param detector the detector
    /// @param context the geometry context
    /// @param show_info boolean to choose if information boxes should be included
    illustrator(const detector_t& detector,
                const geometry_context& context,
                const bool show_info)
        : _detector{detector}, _context{context}, _show_info{show_info} {}

    /// @param detector the detector
    /// @param context the geometry context
    /// @param show_info boolean to choose if information boxes should be included
    /// @param style the styling options to apply
    illustrator(const detector_t& detector,
                const geometry_context& context,
                const bool show_info, const styling::style& style)
        : _detector{detector},
          _context{context},
          _show_info{show_info},
          _style{style} {}

    /// @brief Converts a detray surface in the detector to an svg.
    /// @param identification the id of the svg object.
    /// @param index the index of the surface in the detector.
    /// @param view the display view.
    /// @returns actsvg::svg::object of the detector's surface.
    template <typename view_t>
    inline auto draw_surface(
        const std::string& identification,
        const std::size_t index, const view_t& view) const {
        const auto surface = detray::surface{
            _detector,
            _detector.surface_lookup()[static_cast<detray::dindex>(index)]};
        actsvg::svg::object ret;
        ret._tag = "g";
        ret._id = identification;
        actsvg::svg::object svg_sur;
        std::array<int, 3> color;
        if (surface.is_portal()) {
            auto p_portal = svgtools::conversion::portal<point3_container>(
                _context, _detector, surface);
            svgtools::styling::apply_style(p_portal,
                                           _style._volume_style._portal_style);
            std::copy(p_portal._surface._fill._fc._rgb.begin(),
                      p_portal._surface._fill._fc._rgb.end(), color.begin());
            svg_sur = actsvg::display::portal(identification, p_portal, view);
        } else {
            auto p_surface = svgtools::conversion::surface<point3_container>(
                _context, surface);
            svgtools::styling::apply_style(p_surface,
                                           _style._volume_style._surface_style);
            std::copy(p_surface._fill._fc._rgb.begin(),
                      p_surface._fill._fc._rgb.end(), color.begin());
            svg_sur = actsvg::display::surface(identification, p_surface, view);
        }
        if (_show_info) {
            auto p_information_section =
                svgtools::conversion::information_section<point3>(_context,
                                                                  surface);
            std::copy(color.begin(), color.end(),
                      p_information_section._color.begin());
            ret.add_object(svgtools::meta::display::information_section(
                identification + "_information_section", p_information_section,
                view, _info_screen_offset, svg_sur));
        }
        ret.add_object(svg_sur);
        return ret;
    }

    /// @brief Converts a collection of detray surfaces in the detector to an
    /// svg.
    /// @param identification the id of the svg object.
    /// @param indices the collection of surface indices in the detector to
    /// convert.
    /// @param view the display view.
    /// @returns actsvg::svg::object of the detector's surfaces.
    template <typename iterator_t, typename view_t>
    inline auto draw_surfaces(
        const std::string& identification,
        const iterator_t& indices, const view_t& view) const {
        actsvg::svg::object ret;
        ret._tag = "g";
        ret._id = identification;
        for (const auto index : indices) {
            const auto svg = draw_surface(
                identification + "_surface" + std::to_string(index),
                index, view);
            ret.add_object(svg);
        }
        return ret;
    }

    /// @brief Converts a detray volume in the detector to an svg.
    /// @param identification the id of the svg object.
    /// @param index the index of the volume in the detector.
    /// @param view the display view.
    /// @returns actsvg::svg::object of the detector's volume.
    template <typename view_t>
    inline auto draw_volume(
        const std::string& identification,
        const std::size_t index, const view_t& view) const {
        const auto surface_indices = svgtools::utils::surface_indices(_detector, index);
        return draw_surfaces(identification, surface_indices, view);
    }

    /// @brief Converts a collection of detray volumes in the detector to an
    /// svg.
    /// @param identification the id of the svg object.
    /// @param indices the collection of volume indices in the detector to
    /// convert.
    /// @param view the display view.
    /// @returns actsvg::svg::object of the detector's volumes.
    template <typename iterator_t, typename view_t>
    inline auto draw_volumes(
        const std::string& identification,
        const iterator_t& indices, const view_t& view) const {
        actsvg::svg::object ret;
        ret._tag = "g";
        ret._id = identification;
        for (const auto index : indices) {
            const auto svg =
                draw_volume(identification + "_volume" + std::to_string(index),
                            index, view);
            ret.add_object(svg);
        }
        return ret;
    }

    /// @brief Converts a detray detector to an svg.
    /// @param identification the id of the svg object.
    /// @param view the display view.
    /// @returns actsvg::svg::object of the detector.
    template <typename view_t>
    inline auto draw_detector(
        const std::string& identification,
        const view_t& view) const {
        auto indices =
            detray::views::iota(std::size_t{0u}, _detector.volumes().size());
        return draw_volumes(identification, indices, view);
    }

    /// @brief Converts an intersection record to an svg.
    /// @param identification the id of the svg object.
    /// @param intersection_record the intersection record.
    /// @param view the display view.
    /// @return actsvg::svg::object of the intersection record.
    template <typename view_t>
    inline auto draw_intersections(
        const std::string& identification,
        const std::vector<
            std::pair<detray::dindex,
                      detray::intersection2D<typename detector_t::surface_type,
                                             typename detector_t::transform3>>>&
            intersection_record,
        const view_t& view) const {
        auto p_ir = svgtools::conversion::intersection_record<point3>(
            _context, _detector, intersection_record);
        svgtools::styling::apply_style(p_ir, _style._intersection_style);
        return svgtools::meta::display::intersection_record(identification,
                                                            p_ir, view);
    }

    /// @brief Converts a trajectory to an svg.
    /// @param identification the id of the svg object.
    /// @param trajectory the trajectory (eg. ray or helix).
    /// @param view the display view.
    /// @return actsvg::svg::object of the trajectory.
    template <typename view_t, template <typename> class trajectory_t,
              typename transform3_t>
    inline auto draw_trajectory(const std::string& identification,
                                const trajectory_t<transform3_t>& trajectory,
                                const view_t& view) const {
        auto p_trajectory =
            svgtools::conversion::trajectory<point3>(trajectory);
        svgtools::styling::apply_style(p_trajectory, _style._trajectory_style);
        return svgtools::meta::display::trajectory(identification, p_trajectory,
                                                   view);
    }

    /// @brief Converts a trajectory to an svg.
    /// @param identification the id of the svg object.
    /// @param trajectory the trajectory (eg. ray or helix).
    /// @param view the display view.
    /// @return actsvg::svg::object of the trajectory.
    template <typename view_t, typename point3_container>
    inline auto draw_trajectory(const std::string& identification,
                                const point3_container& points,
                                const view_t& view) const {
        auto p_trajectory =
            svgtools::conversion::trajectory<point3>(points);
        svgtools::styling::apply_style(p_trajectory, _style._trajectory_style);
        return svgtools::meta::display::trajectory(identification, p_trajectory,
                                                   view);
    }

    /// @brief Converts a trajectory and its intersection record to an svg with
    /// a related coloring.
    /// @param identification the id of the svg object.
    /// @param trajectory the trajectory (eg. ray or helix).
    /// @param intersection_record the intersection record.
    /// @param view the display view.
    /// @return actsvg::svg::object of the trajectory and intersection record.
    template <typename view_t, template <typename> class trajectory_t,
              typename transform3_t>
    inline auto draw_intersections_and_trajectory(
        const std::string& identification,
        std::vector<
            std::pair<detray::dindex,
                      detray::intersection2D<typename detector_t::surface_type,
                                             typename detector_t::transform3>>>&
            intersection_record,
        const trajectory_t<transform3_t>& trajectory,
        const view_t& view) const {

        actsvg::svg::object ret;
        ret._tag = "g";
        ret._id = identification;
        auto i_style = svgtools::styling::copy_fill_colors(
            _style._intersection_style, _style._trajectory_style);
        auto p_ir = svgtools::conversion::intersection_record<point3>(
            _context, _detector, intersection_record);
        svgtools::styling::apply_style(p_ir, i_style);
        ret.add_object(svgtools::meta::display::intersection_record(
            identification + "_record", p_ir, view));
        ret.add_object(
            draw_trajectory(identification + "_trajectory", trajectory, view));
        return ret;
    }

    template <typename view_t>
    auto draw_grid(const std::string& identification, const std::size_t index, const view_t& view) const{
        auto p_grid = svgtools::conversion::grid<actsvg::scalar>(_detector, index, view);
        svgtools::styling::apply_style(p_grid, _style._grid_style);
        return actsvg::display::grid(identification, p_grid);
    }


    private:
    using point3 = std::array<actsvg::scalar, 3>;
    using point3_container = std::vector<point3>;
    using geometry_context = typename detector_t::geometry_context;

    const actsvg::point2 _info_screen_offset{-300, 300};
    const detector_t& _detector;
    const geometry_context& _context;
    const bool _show_info = false;
    const styling::style _style = styling::style1;
};

}  // namespace detray::svgtools