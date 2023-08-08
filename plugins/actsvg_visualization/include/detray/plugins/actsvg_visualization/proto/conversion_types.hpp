#pragma once

// Actsvg include(s)
#include "actsvg/core/defs.hpp"
#include "actsvg/proto/detector.hpp"
#include "actsvg/proto/portal.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/proto/volume.hpp"

// System include(s)
#include <array>
#include <vector>

namespace detray::actsvg_visualization::proto {

using point3 = std::array<actsvg::scalar, 3>;
using point3_container = std::vector<point3>;
using proto_surface = actsvg::proto::surface<point3_container>;
using proto_portal = actsvg::proto::portal<point3_container>;
using proto_link = actsvg::proto::portal<point3_container>::link;
using proto_volume = actsvg::proto::volume<point3_container>;
using proto_detector = actsvg::proto::detector<point3_container>;

}  // namespace detray::actsvg_visualization::proto