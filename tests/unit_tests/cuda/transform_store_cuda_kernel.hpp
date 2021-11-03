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

#include "core/transform_store.hpp"

#pragma once

using namespace detray;
using namespace __plugin;

namespace detray {

void transform_test(vecmem::data::vector_view<point3>& input_data,
                    static_transform_store_data& store_data,
                    vecmem::data::vector_view<point3>& output_data);
}
