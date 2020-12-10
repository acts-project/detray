/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "core/surface.hpp"
#include "tools/planar_intersector.hpp"
#include "tests/io/read_csv.hpp"
#include "utils/containers.hpp"

namespace detray
{
    using namespace vector;

    using transform3 = __plugin::transform3;
    using point3 = transform3::point3;
    using vector3 = transform3::vector3;
    using context = transform3::context;

    /** This method creates a number (distances.size()) planes along a direction 
    */
    dvector<surface<transform3>> planes_along_direction(dvector<scalar> distances, vector3 direction)
    {
        context ctx;

        // Rotation matrix
        vector3 z = direction;
        vector3 x = normalize(vector3{0, -z[2], z[1]});

        dvector<surface<transform3>> return_surfaces;
        return_surfaces.reserve(distances.size());
        for (auto &d : distances)
        {
            vector3 t = d * direction;
            transform3 trf(t, z, x, ctx);
            return_surfaces.push_back(surface<transform3>{std::move(trf), std::move(0), std::move(0)});
        }
        return return_surfaces;
    }

} // namespace detray
