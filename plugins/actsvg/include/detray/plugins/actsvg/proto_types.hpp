#pragma once

// Actsvg include(s)
#include "actsvg/proto/surface.hpp"
#include "actsvg/proto/portal.hpp"
#include "actsvg/proto/volume.hpp"
#include "actsvg/proto/detector.hpp"

namespace detray::actsvg_visualization {

using point3 = std::array<actsvg::scalar, 3>;
using point3_container = std::vector<point3>;
using proto_surface = actsvg::proto::surface<point3_container>;
using proto_portal = actsvg::proto::portal<point3_container>;
using proto_link = proto_portal::link;
using proto_volume = actsvg::proto::volume<point3_container>;
using proto_detector = actsvg::proto::detector<point3_container>;


}