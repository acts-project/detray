/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/plugins/svgtools/conversion/point.hpp"
#include "detray/plugins/svgtools/meta/proto/trajectory.hpp"

// System include(s)
#include <vector>

namespace detray::svgtools::conversion {

/// @returns The proto trajectory of a vector of points.
template <typename point3_t>
inline auto trajectory(const std::vector<point3_t>& points) {
    using p_trajectory_t = svgtools::meta::proto::trajectory<point3_t>;
    p_trajectory_t p_trajectory;
    p_trajectory._points = points;
    return p_trajectory;
}
/// @param traj the trajectory
/// @param path_legth the length of the path
/// @param step_size the step size, ie., ds
/// @returns The proto trajectory of a parametrized trajectory
template <typename point3_t, template <typename> class trajectory_t,
          typename transform3_t>
inline auto trajectory(const trajectory_t<transform3_t>& traj, const typename transform3_t::scalar_type path_length = 500.f, const typename transform3_t::scalar_type step_size = 1.f) {
    const typename transform3_t::scalar_type S0 = 0.f;
    std::vector<point3_t> points;
    for (auto s = S0; s < path_length; s += step_size) {
        points.push_back(svgtools::conversion::point<point3_t>(traj.pos(s)));
    }
    return trajectory(points);
}

/// @param traj the trajectory
/// @param path_legth the length of the path
/// @returns The proto trajectory of a ray.
template <typename point3_t, typename transform3_t>
inline auto trajectory(const detray::detail::ray<transform3_t>& traj, const typename transform3_t::scalar_type path_length = 500.f) {
    const typename transform3_t::scalar_type S0 = 0.f;
    std::vector<point3_t> points = {
        svgtools::conversion::point<point3_t>(traj.pos(S0)),
        svgtools::conversion::point<point3_t>(traj.pos(path_length))};
    return trajectory<point3_t>(points);
}

}  // namespace detray::svgtools::conversion