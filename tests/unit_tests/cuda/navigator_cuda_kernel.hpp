/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#if defined(array)
#include "detray/plugins/algebra/array_definitions.hpp"
#elif defined(eigen)
#include "detray/plugins/algebra/eigen_definitions.hpp"
#elif defined(smatrix)
#include "detray/plugins/algebra/smatrix_definitions.hpp"
#elif defined(vc_array)
#include "detray/plugins/algebra/vc_array_definitions.hpp"
#endif

#include "detray/tools/navigator.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"

using namespace detray;

// some useful type declarations
using detector_host_t =
    detector<detector_registry::toy_detector, darray, thrust::tuple,
             vecmem::vector, vecmem::jagged_vector>;
using detector_device_t =
    detector<detector_registry::toy_detector, darray, thrust::tuple,
             vecmem::device_vector, vecmem::jagged_device_vector>;
using navigator_host_t = navigator<detector_host_t>;
using navigator_device_t = navigator<detector_device_t>;
using nav_context = detector_host_t::context;

// detector configuration
constexpr std::size_t n_brl_layers = 4;
constexpr std::size_t n_edc_layers = 3;

// geomery navigation configurations
constexpr unsigned int theta_steps = 100;
constexpr unsigned int phi_steps = 100;

namespace detray {

/// test function for navigator with single state
void navigator_test(
    navigator_view<navigator_host_t> n_data,
    vecmem::data::vector_view<track<nav_context>>& tracks_data,
    vecmem::data::jagged_vector_view<intersection>& candidates_data,
    vecmem::data::jagged_vector_view<dindex>& volume_records_data);

}  // namespace detray
