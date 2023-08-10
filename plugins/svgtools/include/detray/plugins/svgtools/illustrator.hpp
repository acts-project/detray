/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/detector_writer.hpp"
#include "detray/plugins/svgtools/conversion/surface.hpp"
#include "detray/plugins/svgtools/conversion/volume.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"
#include "detray/plugins/svgtools/utils/volume_utils.hpp"
#include "detray/plugins/svgtools/conversion/landmark.hpp"
#include "detray/plugins/svgtools/conversion/intersection_record.hpp"
#include "detray/plugins/svgtools/meta/display/geometry.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"

// System include(s)
#include <array>
#include <vector>

namespace detray::svgtools {

/// @brief SVG generator for a detector and related entities.
template <typename detector_t>
class illustrator {

    public:
    illustrator() = delete;

    illustrator(const detector_t& detector,
                const typename detector_t::name_map& name_map)
        : _detector{detector}, _name_map{name_map} {}

    illustrator(const detector_t& detector,
                const typename detector_t::name_map& name_map,
                const styling::style& style)
        : _detector{detector}, _name_map{name_map}, _style{style} {}

    /// @brief Converts a detray surface in the detector to an svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param index the index of the surface in the detector.
    /// @param view the display view.
    /// @returns actsvg::svg::object of the detector's surface.
    template <typename view_t>
    inline auto draw_surface(const std::string& identification, const typename detector_t::geometry_context& context, const size_t index, const view_t& view) const {
        const auto surface = detray::surface{_detector, _detector.surface_lookup()[static_cast<detray::dindex>(index)]};
        if (surface.is_portal())
        {
            auto p_portal = svgtools::conversion::portal<point3_container>(context, _detector, surface);
            svgtools::styling::apply_style(p_portal, _style._volume_style._portal_style);
            return actsvg::display::portal(identification, p_portal, view);
        }
        auto p_surface = svgtools::conversion::surface<point3_container>(context, surface);
        svgtools::styling::apply_style(p_surface, _style._volume_style._surface_style);
        return actsvg::display::surface(identification, p_surface, view);
    }

    /// @brief Converts a collection of detray surfaces in the detector to an svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param indices the collection of surface indices in the detector to
    /// convert.
    /// @param view the display view.
    /// @returns actsvg::svg::object of the detector's surfaces.
    template <typename iterator_t, typename view_t>
    inline auto draw_surfaces(
        const std::string& identification,
        const typename detector_t::geometry_context& context,
        const iterator_t& indices, const view_t& view) const {
        actsvg::svg::object ret;
        ret._tag = "g";
        ret._id = identification;
        for (const auto index : indices) {
            const auto svg = draw_surface(
                identification + "_surface" + std::to_string(index), context,
                index, view);
            ret.add_object(svg);
        }
        return ret;
    }

    /// @brief Converts a detray volume in the detector to an svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param index the index of the volume in the detector.
    /// @param view the display view.
    /// @returns actsvg::svg::object of the detector's volume.
    template <typename view_t>
    inline auto draw_volume(const std::string& identification, const typename detector_t::geometry_context& context, const size_t index, const view_t& view) const {
        const auto volume = _detector.volume_by_index(static_cast<detray::dindex>(index));
        auto p_volume = svgtools::conversion::volume<point3_container>(context, _detector, volume);
        svgtools::styling::apply_style(p_volume, _style._volume_style);
        return actsvg::display::volume(identification, p_volume, view);
    }

    /// @brief Converts a collection of detray volumes in the detector to an svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param indices the collection of volume indices in the detector to
    /// convert.
    /// @param view the display view.
    /// @returns actsvg::svg::object of the detector's volumes.
    template <typename iterator_t, typename view_t>
    inline auto draw_volumes(
        const std::string& identification,
        const typename detector_t::geometry_context& context,
        const iterator_t& indices, const view_t& view) const {
        actsvg::svg::object ret;
        ret._tag = "g";
        ret._id = identification;
        for (const auto index : indices) {
            const auto svg =
                draw_volume(identification + "_volume" + std::to_string(index),
                            context, index, view);
            ret.add_object(svg);
        }
        return ret;
    }

    /// @brief Converts a detray detector to an svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param detector the detector.
    /// @param view the display view.
    /// @returns actsvg::svg::object of the detector.
    template <typename view_t>
    inline auto draw_detector(
        const std::string& identification,
        const typename detector_t::geometry_context& context,
        const view_t& view) const {
        actsvg::svg::object ret;
        ret._tag = "g";
        ret._id = identification;
        for (size_t index = 0; index < _detector.volumes().size(); index++) {
            const auto svg =
                draw_volume(identification + "_volume" + std::to_string(index),
                            context, index, view);
            ret.add_object(svg);
        }
        return ret;
    }

    /// @brief Converts an intersection record to an svg.
    /// @param identification the id of the svg object.
    /// @param context the geometry context.
    /// @param intersection_record the intersection record.
    /// @param view the display view.
    /// @return actsvg::svg::object of the intersectio record.
    template <typename view_t>
    inline auto draw_intersections(const std::string& identification, const typename detector_t::geometry_context& context, const std::vector<std::pair<detray::dindex, detray::intersection2D<typename detector_t::surface_type, typename detector_t::transform3>>>& intersection_record, const view_t& view) const
    {
        auto p_ir = svgtools::conversion::intersection_record<point3>(context, _detector, intersection_record);
        svgtools::styling::apply_style(p_ir, _style._intersection_style);
        return svgtools::meta::display::intersection_record(identification, p_ir, view);
    }
    
    private:

    private:
    using point3 = std::array<actsvg::scalar, 3>;
    using point3_container = std::vector<point3>;

    const detector_t& _detector;
    const typename detector_t::name_map& _name_map;
    const styling::style _style = styling::style1;
};

}  // namespace detray::svgtools