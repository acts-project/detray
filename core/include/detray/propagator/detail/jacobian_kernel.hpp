/** Detray library, part of the ACTS project
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray::detail {

/// @brief Coordinate frame transformations for track parametrizations
template <template <typename> class frame_t, typename T,
          template <typename> class algebra_t>
struct jacobian_kernel {};

}  // namespace detray::detail
