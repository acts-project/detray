/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project includes(s)
#include "detray/core/detail/single_store.hpp"
#include "detray/definitions/algebra.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>

using namespace detray;
using namespace __plugin;

namespace detray {

using host_transform_store_t =
    single_store<__plugin::transform3<detray::scalar>, vecmem::vector>;

using device_transform_store_t =
    single_store<__plugin::transform3<detray::scalar>, vecmem::device_vector>;

void transform_test(
    vecmem::data::vector_view<point3<detray::scalar> > input_data,
    typename host_transform_store_t::view_type store_data,
    vecmem::data::vector_view<point3<detray::scalar> > output_data,
    std::size_t n_transforms);
}  // namespace detray
