#pragma once

// Project include(s)
#include "detray/plugins/actsvg_visualization/volume.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"

// System include(s)
#include <vector>
#include <iostream>

namespace detray::actsvg_visualization::detector {

struct detector_options{
    volume::volume_options v_options;
};

/// @brief Calculates the proto detector of a detray detector.
///
/// @param d_detector The detray detector.
/// @param context The context.
///
/// @returns An actsvg proto detector representing the detector.
template <typename detector_t>
auto to_proto_detector(
const typename detector_t::geometry_context& context,
const detector_t& d_detector,
const detector_options& d_options
) {
    conversion_types::detector p_detector;
    std::vector<conversion_types::volume> volumes;
    for (const auto& descriptor : d_detector.volumes()){
        std::cout << "doesnt draw middle portals \n";
        auto d_volume = detray::detector_volume{d_detector, descriptor};
        volumes.push_back(volume::to_proto_volume(context, d_detector, d_volume, d_options.v_options));
    }
    p_detector._volumes = volumes;
    return p_detector;
}

template <typename detector_t, typename view_t>
auto to_svg(const typename detector_t::geometry_context& context, const view_t& view, const detector_t& detector, const detector_options& d_options, const std::string& name)
{
    auto p_detector = to_proto_detector(context, detector, d_options);
    return actsvg::display::detector(name, p_detector, view);
}
}