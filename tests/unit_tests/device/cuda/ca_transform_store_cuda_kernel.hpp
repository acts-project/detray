/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project includes(s)
#include "detray/core/detail/context_aware_store.hpp"
#include "detray/definitions/detail/algebra.hpp"

// Vecmem include(s)
#include <vecmem/containers/jagged_device_vector.hpp>

namespace detray {

using algebra_t = ALGEBRA_PLUGIN<scalar>;
using point3 = dpoint3D<algebra_t>;
using transform3 = dtransform3D<algebra_t>;

using host_ca_transform_store_t =
    context_aware_store<transform3, vecmem::jagged_vector,
                        detail::data_context>;

using device_ca_transform_store_t =
    context_aware_store<transform3, vecmem::jagged_device_vector,
                        detail::data_context>;

void ca_transform_test(
    vecmem::data::vector_view<point3> input_data,
    typename host_ca_transform_store_t::view_type ca_store_data,
    vecmem::data::vector_view<point3> output_data, std::size_t n_transforms,
    detail::data_context context);

}  // namespace detray
