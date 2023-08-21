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
inline auto trajectory(const std::vector<point3_t>& points){
    using p_trajectory_t = svgtools::meta::proto::trajectory<point3_t>;
    p_trajectory_t p_trajectory;
    p_trajectory._points = points;
    return p_trajectory;
}

/// @returns The proto trajectory of a vector of a trajectory.
template <typename point3_t, template <typename> class trajectory_t, typename transform3_t>
inline auto trajectory(const trajectory_t<transform3_t>& traj){
    const auto ds = 1.f;
    const auto S0 = 0;
    const auto S1 = 500;
    std::vector<point3_t> points;
    for (float s = S0; s < S1; s+=ds){
        points.push_back(svgtools::conversion::point<point3_t>(traj.pos(s)));
    }
    return trajectory(points);
}

/// @returns The proto trajectory of a vector of a trajectory.
template <typename point3_t, typename transform3_t>
inline auto trajectory(const detray::detail::ray<transform3_t>& traj){
    const auto S0 = 0;
    const auto S1 = 500;
    std::vector<point3_t> points = {svgtools::conversion::point<point3_t>(traj.pos(S0)), svgtools::conversion::point<point3_t>(traj.pos(S1))};
    return trajectory<point3_t>(points);
}

}