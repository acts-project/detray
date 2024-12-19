/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/tutorial/types.hpp"

namespace detray::tutorial {

// Detector
using metadata_t = detray::tutorial::toy_metadata;
using detector_host_t = detector<metadata_t, host_container_types>;
using detector_device_t = detector<metadata_t, device_container_types>;

using mask_id = typename detector_host_t::masks::id;
using acc_id = typename detector_host_t::accel::id;

/// Detector construction tutorial function (prints some detector statistics)
void print(typename detector_host_t::view_type det_data);

}  // namespace detray::tutorial
