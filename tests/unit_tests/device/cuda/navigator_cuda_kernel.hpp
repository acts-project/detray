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

#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>

#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/simulation/track_generators.hpp"

using namespace detray;

using transform3 = __plugin::transform3<scalar>;
using intersection_t = line_plane_intersection;

// some useful type declarations
using detector_host_t = detector<detector_registry::toy_detector, covfie::field,
                                 host_container_types>;
using detector_device_t = detector<detector_registry::toy_detector,
                                   covfie::field_view, device_container_types>;
using navigator_host_t = navigator<detector_host_t>;
using navigator_device_t = navigator<detector_device_t>;
using stepper_t = line_stepper<transform3>;

// detector configuration
constexpr std::size_t n_brl_layers{4u};
constexpr std::size_t n_edc_layers{3u};

// geomery navigation configurations
constexpr unsigned int theta_steps{100u};
constexpr unsigned int phi_steps{100u};

constexpr scalar pos_diff_tolerance{1e-3f};
constexpr scalar overstep_tolerance{-1e-4f};

// dummy propagator state
template <typename navigation_t>
struct prop_state {
    stepper_t::state _stepping;
    navigation_t _navigation;
};

namespace detray {

/// test function for navigator with single state
void navigator_test(
    typename detector_host_t::detector_view_type det_data,
    vecmem::data::vector_view<free_track_parameters<transform3>>& tracks_data,
    vecmem::data::jagged_vector_view<intersection_t>& candidates_data,
    vecmem::data::jagged_vector_view<dindex>& volume_records_data,
    vecmem::data::jagged_vector_view<point3>& position_records_data);

}  // namespace detray
