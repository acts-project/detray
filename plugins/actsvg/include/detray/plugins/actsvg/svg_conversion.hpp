#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg/portal_conversion.hpp"
#include "detray/plugins/actsvg/surface_conversion.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/detector_writer.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization {

    struct surface_style
    {
        const std::vector<actsvg::style::color> fill_colors;
    };

    struct detector_style
    {
        surface_style surface_styling;
        surface_style portal_styling;
    };

    const surface_style default_surface_style{{
        {{0, 128, 0}, 0.5},
        {{34, 139, 34}, 0.5},
        {{60, 179, 113}, 0.5},
        {{85, 107, 47}, 0.5},
        {{124, 252, 0}, 0.5},
        {{154, 205, 50}, 0.5}
        
        }};
    const surface_style default_portal_style{{
        
        {{128, 0, 128}, 0.5},
        {{138, 43, 226}, 0.5},
        {{148, 0, 211}, 0.5},
        {{160, 32, 240}, 0.5},
        {{186, 85, 211}, 0.5},
        {{218, 112, 214}, 0.5}
        
        }};


    const detector_style default_detector_style{default_surface_style, default_portal_style};
    
    template <typename container_t>
    auto pick_random(container_t container){
        int idx = rand() % container.size();
        return container[idx];
    }

    template <typename shape_t, typename view_t>
    actsvg::svg::object svg(
    const std::string& object_name,
    const detray::mask<shape_t>& mask,
    const view_t& view,
    const surface_style& style = default_surface_style
    )
    {
        auto p_surface = convert_mask(mask);
        p_surface._fill = actsvg::style::fill(pick_random(style.fill_colors));
        return actsvg::display::surface(object_name, p_surface, view);
    }

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
            auto fill_color = pick_random(style.portal_styling.fill_colors);
            p_portal._surface._fill = actsvg::style::fill(fill_color);
            return actsvg::display::portal(object_name, p_portal, view);
        }

        auto p_surface = convert_surface(d_surface, context);
        auto fill_color = pick_random(style.surface_styling.fill_colors);
        p_surface._fill = actsvg::style::fill(fill_color);
        return actsvg::display::surface(object_name, p_surface, view);
    }

    template <typename detector_t, typename view_t>
    actsvg::svg::object svg(
    const std::string object_name,
    const detector_t detector,
    const typename detector_t::geometry_context& context,
    const view_t view,
    const detector_style& style = default_detector_style
    )
    {
        actsvg::svg::object ret;
        ret._tag = "g";
        ret._id = object_name;
        for (size_t i = 0; i < detector.surface_lookup().size(); i++){
            const auto description = detector.surface_lookup()[i];
            if (description.volume() != 7){
                continue;
            }
            const auto d_surface = detray::surface{detector, description};
            ret.add_object(svg(object_name + std::to_string(i), detector, d_surface, context, view, style));
        }
        return ret;
    }

    void write_svg(const actsvg::svg::object& svg, const std::string& path){
        actsvg::svg::file file;
        file.add_object(svg);
        detray::io::detail::file_handle stream{path,
                                           "",
                                           std::ios::out | std::ios::trunc};
        *stream << file;
    }

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
    
}  // namespace detray::actsvg_visualization