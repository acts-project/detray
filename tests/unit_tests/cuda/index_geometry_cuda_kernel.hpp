/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#if defined(array)
#include "plugins/algebra/array_definitions.hpp"
#elif defined(eigen)
#include "plugins/algebra/eigen_definitions.hpp"
#elif defined(smatrix)
#include "plugins/algebra/smatrix_definitions.hpp"
#elif defined(vc_array)
#include "plugins/algebra/vc_array_definitions.hpp"
#endif

#include <vecmem/containers/data/jagged_vector_buffer.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

#include "geometry/index_geometry.hpp"

#pragma once

using namespace detray;
using namespace __plugin;

namespace detray {

using geometry = index_ge   ometry<>;
using surface = typename geometry::surface;
using portal = typename geometry::portal;

/// test function for index geometry
void index_geometry_test();

}  // namespace detray
