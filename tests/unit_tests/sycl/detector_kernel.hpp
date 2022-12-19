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
#include "detray/detectors/create_telescope_detector.hpp"
#include "queue_wrapper.hpp"

using namespace detray;

using detector_host_type = detector<detector_registry::telescope_detector,
                                    covfie::field, host_container_types>;

using detector_device_type =
    detector<detector_registry::telescope_detector, covfie::field_view,
             device_container_types>;

namespace detray {

/// test function for propagator with single state
void detector_test(detector_view<detector_host_type> det_data,
                   sycl::queue_wrapper& queue);

}  // namespace detray