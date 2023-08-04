#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/plugins/actsvg_visualization/proto/utils/surface_functors.hpp"


// Actsvg include(s)
#include "actsvg/meta.hpp"
#include "actsvg/core.hpp"

// System include(s)
#include <assert.h>

namespace detray::actsvg_visualization::proto::utils {

/// @brief Checks if the detray surface has a volume link.
template <typename detector_t>
auto has_link(const detray::surface<detector_t>& d_portal){
    const auto d_link_idx = d_portal.template visit_mask<get_link_functor>();
    return d_link_idx != std::numeric_limits<decltype(d_link_idx)>::max();
}
/// @note expects that the detray surface has a volume link.
/// @returns the volume link of the detray surface.
template <typename detector_t>
auto get_link_volume(const detector_t& detector, const detray::surface<detector_t>& d_portal){
    assert(has_link(d_portal));
    const auto d_link_idx = d_portal.template visit_mask<get_link_functor>();
    return detector.volume_by_index(d_link_idx);
}

}