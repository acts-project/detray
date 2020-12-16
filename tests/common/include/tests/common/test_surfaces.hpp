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
            return_surfaces.push_back(surface<transform3>{std::move(trf), 0, 0, false});
        }
        return return_surfaces;
    }

    /** This method creates a barrel description of surfaces */
    dtuple<darray<scalar, 2>, dvector<transform3>> barrel_description(scalar r,
                                                                  scalar stagger_r,
                                                                  unsigned int n_phi,
                                                                  scalar tilt_phi,
                                                                  scalar overlap_rphi,
                                                                  scalar barrel_z,
                                                                  scalar stagger_z,
                                                                  unsigned int n_z)
{
    context ctx;

    // Estimate module dimensions
    scalar module_lx = 2 * r * M_PI * (1 + overlap_rphi) / n_phi;
    scalar module_ly = (barrel_z + (n_z - 1) * stagger_z) / n_z;
    darray<scalar, 2> rectangle_bounds = {0.5 * module_lx, 0.5 * module_ly};

    dvector<transform3> transforms;
    scalar phi_step = 2 * M_PI / n_phi;
    scalar z_step = module_ly - stagger_z;
    scalar z_start = -0.5*n_z*(module_ly - stagger_z);
    for (unsigned int iz = 0; iz < n_z; ++iz)
    {
        scalar z_pos = z_start + iz * z_step;
        for (unsigned int iphi = 0; iphi < n_phi; ++iphi)
        {
            scalar phi = -M_PI + iphi * phi_step;
            scalar r_addon = (iphi%2) ? stagger_r : 0.;
            point3 p = {r_addon * std::cos(phi), r_addon * std::sin(phi), z_pos};
            vector3 z = {std::cos(phi + tilt_phi), std::sin(phi + tilt_phi), 0.};
            vector3 x = {z[1], -z[0], 0.};
            transforms.push_back(transform3(p, z, x, ctx));
        }
    }
    return {rectangle_bounds, transforms};
}

} // namespace detray
