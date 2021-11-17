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

#include "detray/geometry/index_geometry.hpp"

#pragma once

using namespace detray;
using namespace __plugin;

namespace detray {

using geometry = index_geometry<>;
using geometry_device = index_geometry<vecmem::device_vector>;
using surface = typename geometry::surface;
using portal = typename geometry::portal;
using object_id = geometry::object_registry_type::id;

/// test function for index geometry
void index_geometry_test(
    index_geometry_data<geometry>& geometry_data,
    vecmem::data::vector_view<typename geometry::volume_type>& output_data);

}  // namespace detray
