/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/svgtools/conversion/detector.hpp"
#include "detray/plugins/svgtools/conversion/grid.hpp"
#include "detray/plugins/svgtools/conversion/information_section.hpp"
#include "detray/plugins/svgtools/conversion/intersection.hpp"
#include "detray/plugins/svgtools/conversion/landmark.hpp"
#include "detray/plugins/svgtools/conversion/surface.hpp"
#include "detray/plugins/svgtools/conversion/trajectory.hpp"
#include "detray/plugins/svgtools/conversion/volume.hpp"
#include "detray/plugins/svgtools/meta/display/information.hpp"
#include "detray/plugins/svgtools/meta/display/tracking.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"
#include "detray/plugins/svgtools/utils/groups.hpp"
#include "detray/plugins/svgtools/utils/volume_utils.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"

// System include(s)
#include <array>
#include <optional>
#include <vector>

namespace detray::svgtools {

/// @brief SVG generator for a detector and related entities.
///
/// Provides an easy interface for displaying typical objects in a detector.
/// For more flexibility, use the tools in svgtools::conversion to convert
/// detray objects to their respective proto object.
/// Use functions in svgtools::meta::display and actsvg::display to display
/// the proto objects.
///
/// @note Avoid using ids containing spaces or dashes as this seems to cause
/// issues (for instance regarding information boxes). Furthermore, to view
/// information boxes, they must be enabled in the constructor. Furthermore the
/// svg viewer (opening the file after it is created) must support animations.
template <typename detector_t>
class illustrator {

    using point3 = std::array<actsvg::scalar, 3>;
    using point3_container = std::vector<point3>;

    public:
    illustrator() = delete;

    /// @param detector the detector
    /// @param name_map naming scheme of the detector
    /// @param style display style
    ///
    /// @note information boxes are enabled by default
    illustrator(const detector_t& detector,
                const typename detector_t::name_map& name_map,
                const styling::style& style =
                    detray::svgtools::styling::tableau_colorblind::style)
        : _detector{detector}, _name_map{name_map}, _style{style} {}

    /// @returns the detector and volume names
    const typename detector_t::name_map& names() { return _name_map; }

    /// Access the illustrator style
    styling::style& style() { return _style; }

    /// Toggle info boxes
    void show_info(bool toggle = true) { _show_info = toggle; }
    /// Toggle info boxes
    void hide_grids(bool toggle = true) { _hide_grids = toggle; }
    /// Toggle portal surfaces
    void hide_portals(bool toggle = true) { _hide_portals = toggle; }
    /// Toggle passive surfaces
    void hide_passives(bool toggle = true) { _hide_passives = toggle; }

    /// @brief Converts a single detray surface of the detector to an svg.
    ///
    /// @param index the index of the surface in the detector.
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @returns @c actsvg::svg::object of the detector's surface.
    template <typename view_t>
    inline auto draw_surface(
        const std::size_t index, const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {

        const auto surface =
            detray::surface{_detector, static_cast<detray::dindex>(index)};

        actsvg::svg::object ret;
        const auto& style = _style._detector_style._volume_style;

        if (surface.is_portal()) {
            auto p_portal = svgtools::conversion::portal<point3_container>(
                gctx, _detector, surface, false);

            svgtools::styling::apply_style(p_portal, style._portal_style);

            // Draw the portal directly
            std::string id = p_portal._name + "_" + svg_id(view);
            ret = actsvg::display::portal(std::move(id), p_portal, view);
        } else {
            auto p_surface =
                svgtools::conversion::surface<point3_container>(gctx, surface);

            svgtools::styling::apply_style(p_surface, style._surface_style);

            // Draw the surface directly
            std::string id = p_surface._name + "_" + svg_id(view);
            ret = actsvg::display::surface(std::move(id), p_surface, view);
        }
        // Add an optional info box
        if (_show_info) {
            auto p_information_section =
                svgtools::conversion::information_section<point3>(gctx,
                                                                  surface);

            auto info_box = svgtools::meta::display::information_section(
                ret._id + "_info_box", p_information_section, view,
                _info_screen_offset, ret);
            ret.add_object(info_box);
        }

        return ret;
    }

    /// @brief Converts a multiple of detray surfaces of the detector to an svg.
    ///
    /// @param indices the collection of surface indices in the detector to
    /// convert.
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @returns @c actsvg::svg::object of the detector's surfaces.
    template <typename iterator_t, typename view_t>
    inline auto draw_surfaces(
        const iterator_t& indices, const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {

        auto ret = svgtools::utils::group(_name_map.at(0) + "_surfaces_" +
                                          svg_id(view));

        for (const auto index : indices) {
            const auto svg = draw_surface(index, view, gctx);
            ret.add_object(svg);
        }

        return ret;
    }

    /// @brief Converts a detray volume of the detector to an svg.
    ///
    /// @param index the index of the volume in the detector.
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @returns @c actsvg::svg::object of the detector's volume.
    template <typename view_t>
    inline auto draw_volume(
        const std::size_t index, const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {

        const auto d_volume =
            detector_volume{_detector, static_cast<detray::dindex>(index)};

        auto p_volume = svgtools::conversion::volume<point3_container>(
            gctx, _detector, d_volume, view, _hide_portals, _hide_passives,
            _hide_grids);

        // Apply the specific styling
        svgtools::styling::apply_style(p_volume,
                                       _style._detector_style._volume_style);

        std::string id = d_volume.name(_name_map) + "_" + svg_id(view);
        auto vol_svg = actsvg::display::volume(std::move(id), p_volume, view);

        if (!_hide_grids) {
            auto grid_svg =
                actsvg::display::grid(id + "_grid", p_volume._surface_grid);
            vol_svg.add_object(grid_svg);
        }

        return vol_svg;
    }

    /// @brief Converts multiple detray volumes of the detector to an svg.
    ///
    /// @param indices the collection of volume indices in the detector to
    /// convert.
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @returns @c actsvg::svg::object of the detector's volumes.
    template <typename range_t, typename view_t>
    inline auto draw_volumes(
        const range_t& indices, const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {

        auto ret = svgtools::utils::group(_name_map.at(0) + "_volumes_" +
                                          svg_id(view));
        for (const auto index : indices) {
            ret.add_object(draw_volume(index, view, gctx));
        }

        return ret;
    }

    /// @brief Converts a detray detector to an svg.
    ///
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @returns @c actsvg::svg::object of the detector.
    template <typename view_t>
    inline auto draw_detector(
        const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {

        auto p_detector = svgtools::conversion::detector<point3_container>(
            gctx, _detector, view, _hide_portals, _hide_passives, _hide_grids);

        // Apply the specific styling
        svgtools::styling::apply_style(p_detector, _style._detector_style);

        std::string id = _name_map.at(0) + "_" + svg_id(view);

        return actsvg::display::detector(std::move(id), p_detector, view);
    }

    /// @brief Converts a point to an svg.
    ///
    /// @param prefix the id of the svg object.
    /// @param point the point.
    /// @param view the display view.
    ///
    /// @return actsvg::svg::object of the point.
    template <typename view_t, typename point_t>
    inline auto draw_landmark(const std::string& prefix, const point_t& point,
                              const view_t& view) const {
        auto p_landmark = svgtools::conversion::landmark<point3>(point);
        svgtools::styling::apply_style(p_landmark, _style._landmark_style);
        return svgtools::meta::display::landmark(prefix + "_landmark",
                                                 p_landmark, view);
    }

    /// @brief Converts an intersection record to an svg.
    ///
    /// @param prefix the id of the svg object.
    /// @param intersection_record the intersection record.
    /// @param dir the direction of the trajectory.
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @return @c actsvg::svg::object of the intersectio record.
    template <typename view_t, typename intersection_t>
    inline auto draw_intersections(
        const std::string& prefix,
        const std::vector<std::pair<detray::dindex, intersection_t>>&
            intersection_record,
        const typename detector_t::vector3 dir, const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {

        dvector<intersection_t> intersections{};
        intersections.reserve(intersection_record.size());

        for (auto& ir : intersection_record) {
            intersections.push_back(ir.second);
        }

        return draw_intersections(prefix, intersections, dir, view, gctx);
    }

    /// @brief Converts a collection of intersections to an svg.
    ///
    /// @param prefix the id of the svg object.
    /// @param intersections the intersection collection.
    /// @param dir the direction of the trajectory.
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @return @c actsvg::svg::object of the intersectio record.
    template <typename view_t, typename intersection_t>
    inline auto draw_intersections(
        const std::string& prefix, const dvector<intersection_t>& intersections,
        const typename detector_t::vector3 dir, const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {

        auto p_ir = svgtools::conversion::intersection<point3>(
            _detector, intersections, dir, gctx);

        svgtools::styling::apply_style(p_ir, _style._intersection_style);

        return svgtools::meta::display::intersection(prefix, p_ir, view);
    }

    /// @brief Converts a trajectory to an svg.
    ///
    /// @param prefix the id of the svg object.
    /// @param trajectory the trajectory (eg. ray or helix).
    /// @param view the display view.
    ///
    /// @return @c actsvg::svg::object of the trajectory.
    template <typename view_t, template <typename> class trajectory_t,
              typename transform3_t>
    inline auto draw_trajectory(const std::string& prefix,
                                const trajectory_t<transform3_t>& trajectory,
                                const typename transform3_t::scalar_type path,
                                const view_t& view) const {
        auto p_trajectory =
            svgtools::conversion::trajectory<point3>(trajectory, path);

        svgtools::styling::apply_style(p_trajectory, _style._trajectory_style);

        return svgtools::meta::display::trajectory(prefix + "_traj",
                                                   p_trajectory, view);
    }

    /// @brief Converts a trajectory to an svg.
    ///
    /// @param prefix the id of the svg object.
    /// @param trajectory the trajectory (eg. ray or helix).
    /// @param view the display view.
    ///
    /// @return actsvg::svg::object of the trajectory.
    template <typename view_t, typename point3_container>
    inline auto draw_trajectory(const std::string& prefix,
                                const point3_container& points,
                                const view_t& view) const {

        auto p_trajectory = svgtools::conversion::trajectory<point3>(points);

        svgtools::styling::apply_style(p_trajectory, _style._trajectory_style);

        return svgtools::meta::display::trajectory(prefix + "_traj",
                                                   p_trajectory, view);
    }

    /// @brief Converts a trajectory and its intersection record to an svg with
    /// a related coloring.
    ///
    /// @param prefix the id of the svg object.
    /// @param trajectory the trajectory (eg. ray or helix).
    /// @param intersection_record the intersection record.
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @return @c actsvg::svg::object of the trajectory and intersection
    /// record.
    template <typename view_t, class trajectory_t, class intersection_t>
    inline auto draw_intersections_and_trajectory(
        const std::string& prefix,
        const std::vector<std::pair<detray::dindex, intersection_t>>&
            intersection_record,
        const trajectory_t& trajectory, const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {

        dvector<intersection_t> intersections{};
        intersections.reserve(intersection_record.size());
        for (auto& ir : intersection_record) {
            intersections.push_back(ir.second);
        }

        return draw_intersections_and_trajectory(prefix, intersections,
                                                 trajectory, view, gctx);
    }

    /// @brief Converts a trajectory and its intersections to an svg with a
    /// related coloring.
    ///
    /// @param prefix the id of the svg object.
    /// @param intersections the intersection record.
    /// @param trajectory the trajectory (eg. ray or helix).
    /// @param view the display view.
    /// @param max_path the maximal path length from trajectory origin.
    /// @param gctx the geometry context.
    ///
    /// @return @c actsvg::svg::object of the trajectory and intersections
    template <typename view_t, class trajectory_t, typename intersection_t>
    inline auto draw_intersections_and_trajectory(
        const std::string& prefix, const dvector<intersection_t>& intersections,
        const trajectory_t& trajectory, const view_t& view,
        typename detector_t::scalar_type max_path = 500.f,
        const typename detector_t::geometry_context& gctx = {}) const {

        actsvg::svg::object ret;
        ret._tag = "g";
        ret._id = prefix;

        if (not intersections.empty()) {

            auto p_ir = svgtools::conversion::intersection<point3>(
                _detector, intersections, trajectory.dir(0.f), gctx);
            svgtools::styling::apply_style(p_ir, _style._landmark_style);

            ret.add_object(svgtools::meta::display::intersection(
                prefix + "_record", p_ir, view));

            // The intersection record is always sorted by path length
            const auto sf{
                detray::surface{_detector, intersections.back().sf_desc}};
            const auto final_pos = sf.local_to_global(
                gctx, intersections.back().local, trajectory.dir(0.f));
            max_path = getter::norm(final_pos - trajectory.pos(0.f));
        }

        ret.add_object(draw_trajectory(prefix + "_trajectory", trajectory,
                                       max_path, view));
        return ret;
    }

    private:
    /// @returns the string id of a view
    template <typename view_t>
    std::string svg_id(const view_t& view) const {
        return view._axis_names[0] + view._axis_names[1];
    }

    const actsvg::point2 _info_screen_offset{-300, 300};
    const detector_t& _detector;
    const typename detector_t::name_map& _name_map;
    bool _show_info = true;
    bool _hide_grids = false;
    bool _hide_portals = false;
    bool _hide_passives = false;
    styling::style _style = styling::tableau_colorblind::style;
};

}  // namespace detray::svgtools
