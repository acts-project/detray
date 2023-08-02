#pragma once

// Actsvg include(s)
#include "actsvg/proto/surface.hpp"
#include "actsvg/proto/portal.hpp"
#include "actsvg/proto/volume.hpp"
#include "actsvg/proto/detector.hpp"

namespace detray::actsvg_visualization::conversion_types {

using point3 = std::array<actsvg::scalar, 3>;
using point3_container = std::vector<point3>;
using surface = actsvg::proto::surface<point3_container>;
using portal = actsvg::proto::portal<point3_container>;
using link = actsvg::proto::portal<point3_container>::link;
using volume = actsvg::proto::volume<point3_container>;
using detector = actsvg::proto::detector<point3_container>;

}