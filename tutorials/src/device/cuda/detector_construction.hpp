/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/tutorial/types.hpp"

// Covfie include(s).
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>

namespace detray::tutorial {

// Detector
using detector_host_t = detector<detector_registry::toy_detector,
                                 covfie::field<const_bfield_bknd_t>, host_container_types>;
using detector_device_t = detector<detector_registry::toy_detector,
                                   covfie::field_view<const_bfield_bknd_t>, device_container_types>;

using mask_id = typename detector_host_t::masks::id;

/// Detector construction tutorial function (prints some detector statistics)
void print(typename detector_host_t::detector_view_type det_data);

}  // namespace detray::tutorial
