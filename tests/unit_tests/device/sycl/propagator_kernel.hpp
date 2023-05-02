/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "queue_wrapper.hpp"
#include "tests/common/test_base/propagator_test.hpp"

namespace detray {

/// test function for propagator with single state
void propagator_test(
    typename detector_host_type::detector_view_type det_data,
    vecmem::data::vector_view<free_track_parameters<transform3>> &tracks_data,
    vecmem::data::jagged_vector_view<intersection_t> &candidates_data,
    vecmem::data::jagged_vector_view<scalar> &path_lengths_data,
    vecmem::data::jagged_vector_view<vector3> &positions_data,
    vecmem::data::jagged_vector_view<free_matrix> &jac_transports_data,
    sycl::queue_wrapper queue);

}  // namespace detray
