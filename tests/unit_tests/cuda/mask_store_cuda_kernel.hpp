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

#include "core/mask_store.hpp"
#include "masks/masks.hpp"
#include "vecmem/utils/cuda/copy.hpp"

#pragma once

using namespace detray;
using namespace __plugin;

const int n_points = 1000;

namespace detray {

using annulus = annulus2<>;
using cylinder = cylinder3<>;
using rectangle = rectangle2<>;
using ring = ring2<>;
using single = single3<0>;
using trapezoid = trapezoid2<>;

void mask_test(mask_store_data<rectangle, trapezoid, ring, cylinder, single,
                               annulus>& store_data,
               vecmem::data::vector_view<point2>& input_point2_data,
               vecmem::data::vector_view<point3>& input_point3_data,
               vecmem::data::jagged_vector_view<int>& output_data);

}  // namespace detray
