/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <climits>
#include <optional>
#include <tuple>

namespace detray
{

    /** Track struct for navigation through the detector */
    template <typename transform_type>
    struct track
    {
        using point3 = typename transform_type::point3;
        using vector3 = typename transform_type::vector3;
        using context = typename transform_type::context;

        point3 pos = {0., 0., 0.};
        vector3 dir = {0., 0., 0.};
        vector3 bfield = {0., 0., 0.};
        scalar momentum = std::numeric_limits<scalar>::infinity();
        context ctx;

        scalar overstep_tolerance = 0.;
    };

} // namespace detray
