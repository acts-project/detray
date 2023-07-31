#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/plugins/actsvg/volume_conversion.hpp"

// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/proto/detector.hpp"

// System include(s)
#include <type_traits>
#include <vector>

namespace detray::actsvg_visualization {

/// @brief Calculates the proto detector of a detray detector.
///
/// @param d_detector The detray detector.
/// @param context The context.
///
/// @returns An actsvg proto detector representing the detector.
template <typename detector_t>
auto convert_detector(
    const detector_t& d_detector,
    const typename detector_t::geometry_context& context) {
        proto_detector p_detector;
        std::vector<proto_volume> volumes;
        for (const auto& descriptor : d_detector.volumes()){
            auto d_volume = detray::detector_volume{d_detector, descriptor};
            volumes.push_back(convert_volume(d_detector, d_volume, context));
        }
        p_detector._volumes = volumes;
        return p_detector;
}
}