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
#include "detray/definitions/indexing.hpp"
#include "detray/utils/ranges.hpp"

// vecmem core
#include "vecmem/containers/device_vector.hpp"

namespace detray {

/// Test @c detray::views::single
void single(const dindex value, dindex& check);

/*void iota(vecmem::data::vector_view<dindex>& check_data,
          const darray<dindex, 2> interval);

struct uint_holder {
    unsigned int ui = 0;
};

void enumerate(,
                vecmem::data::vector_view<dindex>& check_data);

void chain(vecmem::data::vector_view<dindex>& idx_data,
                        vecmem::data::vector_view<unsigned int>& uint_data,
                        vecmem::data::vector_view<uint_holder>& seq_data);

void pick(vecmem::data::vector_view<int>& check_data,
                   vecmem::data::vector_view<int>& seq_data,
                   const size_t& begin, const size_t& end);

void subrange(vecmem::data::vector_view<int>& check_data,
                   vecmem::data::vector_view<int>& seq_data,
                   const size_t& begin, const size_t& end);*/

}  // namespace detray
