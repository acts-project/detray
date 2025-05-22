/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/propagation_config.hpp"

// Detray test include(s)
#include "detray/test/common/track_generators.hpp"
#include "detray/test/framework/types.hpp"

namespace detray {

// some useful type declarations
using metadata_t = test::toy_metadata;
using test_algebra = metadata_t::algebra_type;
using scalar = dscalar<test_algebra>;
using point3 = dpoint3D<test_algebra>;
using detector_host_t = detector<metadata_t, host_container_types>;
using detector_device_t = detector<metadata_t, device_container_types>;

using intersection_t =
    intersection2D<typename detector_device_t::surface_type, test_algebra>;

using navigator_host_t = navigator<detector_host_t>;
using navigator_device_t = navigator<detector_device_t>;
using stepper_t = line_stepper<test_algebra>;

// detector configuration
constexpr std::size_t n_brl_layers{4u};
constexpr std::size_t n_edc_layers{3u};

// geomery navigation configurations
constexpr unsigned int theta_steps{100u};
constexpr unsigned int phi_steps{100u};

constexpr scalar pos_diff_tolerance{1e-3f};

// dummy propagator state
template <typename navigation_t>
struct prop_state {
    using context_t = typename navigation_t::detector_type::geometry_context;
    stepper_t::state _stepping;
    navigation_t _navigation;
    context_t _context{};
};

/// test function for navigator with single state
void navigator_test(
    typename detector_host_t::view_type det_data, navigation::config& nav_cfg,
    stepping::config& step_cfg,
    vecmem::data::vector_view<free_track_parameters<test_algebra>>& tracks_data,
    vecmem::data::jagged_vector_view<dindex>& volume_records_data,
    vecmem::data::jagged_vector_view<point3>& position_records_data);

}  // namespace detray
