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
#include "detray/utils/ranges.hpp"

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
                    detray::svgtools::styling::tableau_colorblind)
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

    /// @brief Converts a detray surface in the detector to an svg.
    ///
    /// @param prefix the id of the svg object.
    /// @param index the index of the surface in the detector.
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @returns @c actsvg::svg::object of the detector's surface.
    template <typename view_t>
    inline auto draw_surface(
        const std::string& prefix, const std::size_t index, const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {

        const auto surface = detray::surface{
            _detector,
            _detector.surface_lookup()[static_cast<detray::dindex>(index)]};

        auto ret = svgtools::utils::group(prefix);
        actsvg::svg::object svg_sur;
        std::array<int, 3> color;

        if (surface.is_portal()) {
            if (!_hide_portals) {
                auto p_portal = svgtools::conversion::portal<point3_container>(
                    gctx, _detector, surface);

                svgtools::styling::apply_style(
                    p_portal, _style._volume_style._portal_style,
                    _style._do_random_coloring);

                std::copy(p_portal._surface._fill._fc._rgb.begin(),
                          p_portal._surface._fill._fc._rgb.end(),
                          color.begin());

                svg_sur = actsvg::display::portal(prefix, p_portal, view);
            }
        } else if (!(surface.is_passive() && _hide_passives)) {
            auto p_surface =
                svgtools::conversion::surface<point3_container>(gctx, surface);

            svgtools::styling::apply_style(p_surface,
                                           _style._volume_style._surface_style,
                                           _style._do_random_coloring);

            std::copy(p_surface._fill._fc._rgb.begin(),
                      p_surface._fill._fc._rgb.end(), color.begin());

            svg_sur = actsvg::display::surface(prefix, p_surface, view);
        }
        if (_show_info &&
            !(_hide_portals && _hide_passives && !surface.is_sensitive())) {
            auto p_information_section =
                svgtools::conversion::information_section<point3>(gctx,
                                                                  surface);
            std::copy(color.begin(), color.end(),
                      p_information_section._color.begin());

            ret.add_object(svgtools::meta::display::information_section(
                prefix + "_sf_" + std::to_string(index), p_information_section,
                view, _info_screen_offset, svg_sur));
        }

        ret.add_object(svg_sur);
        return ret;
    }

    /// @brief Converts a collection of detray surfaces in the detector to an
    /// svg.
    ///
    /// @param prefix the id of the svg object.
    /// @param indices the collection of surface indices in the detector to
    /// convert.
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @returns @c actsvg::svg::object of the detector's surfaces.
    template <typename iterator_t, typename view_t>
    inline auto draw_surfaces(
        const std::string& prefix, const iterator_t& indices,
        const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {
        auto ret = svgtools::utils::group(prefix);
        for (const auto index : indices) {
            const auto svg = draw_surface(prefix, index, view, gctx);
            ret.add_object(svg);
        }
        return ret;
    }

    /// @brief Converts a surface grid to svg.
    ///
    /// @param prefix the id of the svg object.
    /// @param index the index of the grid's volume
    /// @param view the display view.
    ///
    /// @return @c actsvg::svg::object of the trajectory and intersections
    template <typename view_t>
    std::optional<actsvg::svg::object> draw_grid(const std::string& prefix,
                                                 const std::size_t index,
                                                 const view_t& view) const {

        if (auto p_grid_ptr = svgtools::conversion::grid<actsvg::scalar>(
                _detector, index, view)) {

            svgtools::styling::apply_style(*p_grid_ptr, _style._grid_style);

            std::string name{prefix + "_" +
                             _name_map.at(static_cast<detray::dindex>(index)) +
                             "_grid"};
            return actsvg::display::grid(std::move(name), *p_grid_ptr);
        }

        return {};
    }

    /// @brief Converts a detray volume in the detector to an svg.
    ///
    /// @param prefix the id of the svg object.
    /// @param index the index of the volume in the detector.
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @returns @c actsvg::svg::object of the detector's volume.
    template <typename view_t>
    inline auto draw_volume(
        const std::string& prefix, const std::size_t index, const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {

        auto ret = svgtools::utils::group(
            prefix + "_" + _name_map.at(static_cast<detray::dindex>(index)));

        const auto surface_indices =
            svgtools::utils::surface_indices(_detector, index);
        const auto sf_svg = draw_surfaces(prefix, surface_indices, view, gctx);

        ret.add_object(sf_svg);

        if (!_hide_grids) {
            const auto grid_svg = draw_grid(prefix, index, view);
            if (grid_svg.has_value()) {
                ret.add_object(grid_svg.value());
            }
        }

        return ret;
    }

    /// @brief Converts a collection of detray volumes in the detector to an
    /// svg.
    ///
    /// @param prefix the id of the svg object.
    /// @param indices the collection of volume indices in the detector to
    /// convert.
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @returns @c actsvg::svg::object of the detector's volumes.
    template <typename range_t, typename view_t>
    inline auto draw_volumes(
        const std::string& prefix, const range_t& indices, const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {

        auto ret = svgtools::utils::group(prefix);
        for (const auto index : indices) {
            const auto vol_svg = draw_volume(prefix, index, view, gctx);
            ret.add_object(vol_svg);
        }

        return ret;
    }

    /// @brief Converts a detray detector to an svg.
    ///
    /// @param prefix the id of the svg object.
    /// @param view the display view.
    /// @param gctx the geometry context.
    ///
    /// @returns @c actsvg::svg::object of the detector.
    template <typename view_t>
    inline auto draw_detector(
        const std::string& prefix, const view_t& view,
        const typename detector_t::geometry_context& gctx = {}) const {
        auto indices =
            detray::views::iota(std::size_t{0u}, _detector.volumes().size());
        return draw_volumes(prefix, indices, view, gctx);
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
        svgtools::styling::apply_style(p_landmark, _style._landmark_style,
                                       _style._do_random_coloring);
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

        svgtools::styling::apply_style(p_ir, _style._intersection_style,
                                       _style._do_random_coloring);

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

        svgtools::styling::apply_style(p_trajectory, _style._trajectory_style,
                                       _style._do_random_coloring);

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
        svgtools::styling::apply_style(p_trajectory, _style._trajectory_style,
                                       _style._do_random_coloring);
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
            auto i_style = svgtools::styling::copy_fill_colors(
                _style._intersection_style, _style._trajectory_style);
            i_style._fill_colors[0]._opacity = 0.5f;
            i_style._marker_type = "o";

            auto p_ir = svgtools::conversion::intersection<point3>(
                _detector, intersections, trajectory.dir(0.f), gctx);
            svgtools::styling::apply_style(p_ir, i_style,
                                           _style._do_random_coloring);

            ret.add_object(svgtools::meta::display::intersection(
                prefix + "_record", p_ir, view));

            // The intersection record is always sorted by path length
            // const auto max_path{intersections.back().path};
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
    const actsvg::point2 _info_screen_offset{-300, 300};
    const detector_t& _detector;
    const typename detector_t::name_map& _name_map;
    bool _show_info = true;
    bool _hide_grids = false;
    bool _hide_portals = false;
    bool _hide_passives = false;
    styling::style _style = styling::style1;
};

}  // namespace detray::svgtools
