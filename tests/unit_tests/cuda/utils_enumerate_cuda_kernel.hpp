/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// detray test
#include "tests/common/test_defs.hpp"

// detray core
#include "detray/utils/enumerate.hpp"
#include "detray/utils/indexing.hpp"

// vecmem core
#include "vecmem/containers/device_vector.hpp"

namespace detray {

// test function for enumeration with single integer
void sequence_single(vecmem::data::vector_view<dindex>& check_data,
                     vecmem::data::vector_view<dindex>& single_data);

// test function for enumeration with range
void sequence_range(const darray<dindex, 2> range,
                    vecmem::data::vector_view<dindex>& check_data);

}  // namespace detray
