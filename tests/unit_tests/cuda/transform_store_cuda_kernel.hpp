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

#include "detray/core/transform_store.hpp"

#pragma once

using namespace detray;
using namespace __plugin;

namespace detray {

void transform_test(
    vecmem::data::vector_view<point3<detray::scalar> >& input_data,
    static_transform_store_data<static_transform_store<> >& store_data,
    vecmem::data::vector_view<point3<detray::scalar> >& output_data);
}
