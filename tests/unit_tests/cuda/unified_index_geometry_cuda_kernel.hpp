
/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#if defined(array)
#include "detray/plugins/algebra/array_definitions.hpp"
#elif defined(eigen)
#include "detray/plugins/algebra/eigen_definitions.hpp"
#elif defined(smatrix)
#include "detray/plugins/algebra/smatrix_definitions.hpp"
#elif defined(vc_array)
#include "detray/plugins/algebra/vc_array_definitions.hpp"
#endif

#include <vecmem/containers/data/jagged_vector_buffer.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

#include "detray/geometry/unified_index_geometry.hpp"

#pragma once

using namespace detray;
using namespace __plugin;

namespace detray {

using geometry = unified_index_geometry<>;
using surface = typename geometry::surface;
using portal = typename geometry::portal;

/// test function for index geometry
void unified_index_geometry_test(
    unified_index_geometry_data<geometry>& geometry_data,
    vecmem::data::vector_view<typename geometry::volume_type>& output_data);

}  // namespace detray
