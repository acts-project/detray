/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <string>

namespace detray {

/// @brief Wraps any shape, but does not enforce boundaries
template <typename shape_t>
class unbounded {
    public:
    using shape = shape_t;
    using boundaries = typename shape::boundaries;

    /// The name for this shape
    inline static const std::string name = "unbounded " + shape::name;

    /// The measurement dimension
    inline static constexpr const std::size_t meas_dim = shape::meas_dim;

    /// Local coordinate frame for boundary checks
    template <typename algebra_t>
    using local_frame_type =
        typename shape::template local_frame_type<algebra_t>;
    /// Local point type (2D)
    template <typename algebra_t>
    using loc_point_type = typename shape::template loc_point_type<algebra_t>;

    /// Measurement frame
    template <typename algebra_t>
    using measurement_frame_type =
        typename shape::template measurement_frame_type<algebra_t>;
    /// Local measurement point (2D)
    template <typename algebra_t>
    using measurement_point_type =
        typename shape::template measurement_point_type<algebra_t>;

    /// Underlying surface geometry
    template <typename algebra_t>
    using intersector_type =
        typename shape::template intersector_type<algebra_t>;

    /// @brief Check boundary values for a local point.
    ///
    /// @tparam bounds_t any type of boundary values
    ///
    /// @note the parameters are ignored
    ///
    /// @return always true
    template <typename bounds_t, typename point_t, typename scalar_t>
    DETRAY_HOST_DEVICE inline constexpr bool check_boundaries(
        const bounds_t& /*bounds*/, const point_t& /*loc_p*/,
        const scalar_t /*tol*/) const {
        return true;
    }

    template <typename param_t>
    DETRAY_HOST_DEVICE inline typename param_t::point2 to_measurement(
        param_t& param, const typename param_t::point2& offset = {0, 0}) const {
        return param.local() + offset;
    }
};

}  // namespace detray