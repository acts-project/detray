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

    /** This method creates a number (distances.size()) planes along a direction 
    */
    dvector<surface<transform3>> planes_along_direction(dvector<scalar> distances, vector3 direction)
    {
        // Rotation matrix
        vector3 z = direction;
        vector3 x = normalize(vector3{0, -z[2], z[1]});

        dvector<surface<transform3>> return_surfaces;
        return_surfaces.reserve(distances.size());
        for (auto &d : distances)
        {
            vector3 t = d * direction;
            transform3 trf(t, z, x);
            return_surfaces.push_back(surface<transform3>{std::move(trf), 0, 0, false});
        }
        return return_surfaces;
    }


    /** This method creates a barrel description of surfaces
     * 
     * @param r_inner The inner radius of the endcap disc
     * @param r_outer The outer radius of the endcap disc
     * @param z_pos The nominal z position of the endcap disc
     * @param stagger_z The possible staggering in z 
     * @param n_phi_half The half number of modules in phi
     * @param overlap_phi The overlap in phi coverage (approximated)
     * 
     * @returns a tuyple for rectangle descriptions and transforms
     */
    dtuple<darray<scalar, 3>, dvector<transform3>> endcap_description(scalar r_inner,
                                                                      scalar r_outer,
                                                                      scalar z_pos,
                                                                      scalar stagger_z,
                                                                      unsigned int n_phi_half,
                                                                      scalar overlap_phi)
    {
        scalar module_inner_lx = r_inner*M_PI *(1+overlap_phi)/(n_phi_half);
        scalar module_outer_lx = r_inner*M_PI *(1+overlap_phi)/(n_phi_half);
        scalar module_hy = 0.5*(r_outer - r_inner);

        darray<scalar, 3> trapezoid_values = { 0.5*module_inner_lx, 0.5*module_outer_lx, module_hy };
        scalar phi_step = M_PI/n_phi_half;
        dvector<transform3> transforms;

        scalar r = 0.5 * (r_inner + r_outer );

        for (unsigned int iphi = 0; iphi < 2 * n_phi_half;  ++iphi){
            scalar phi = -M_PI + iphi * phi_step;
            scalar z_addon = (iphi % 2) ? -stagger_z : stagger_z;
            scalar cos_phi = std::cos(phi);
            scalar sin_phi = std::sin(phi);
            point3 p = {r * cos_phi, r * sin_phi, z_pos + z_addon};           
            vector3 z = {0., 0., 1.};
            vector3 x = {-sin_phi, cos_phi, 0.};
            transforms.push_back(transform3(p, z, x));
        }
        return { trapezoid_values, transforms };

    }

    /** This method creates a barrel description of surfaces
     * 
     * @param r The radius of the barrel cylinder
     * @param stagger_r The (optional) staggering in r
     * @param n_phi The number of modules in phi
     * @param tilt_phi The tilt of the modules in phi
     * @param overlap_rphi The overlap in rphi coverage (approximated)
     * @param barrel_z The length of the barrel in z
     * @param stagger_z The (optional) staggering in z
     * 
     * @returns a tuyple for rectangle descriptions and transforms
     */
    dtuple<darray<scalar, 2>, dvector<transform3>> barrel_description(scalar r,
                                                                      scalar stagger_r,
                                                                      unsigned int n_phi,
                                                                      scalar tilt_phi,
                                                                      scalar overlap_rphi,
                                                                      scalar barrel_z,
                                                                      scalar stagger_z,
                                                                      unsigned int n_z)
    {
        // Estimate module dimensions
        scalar module_lx = 2 * r * M_PI * (1 + overlap_rphi) / n_phi;
        scalar module_ly = (barrel_z + (n_z - 1) * stagger_z) / n_z;
        darray<scalar, 2> rectangle_bounds = {0.5 * module_lx, 0.5 * module_ly};

        dvector<transform3> transforms;
        scalar phi_step = 2 * M_PI / n_phi;
        scalar z_step = module_ly - stagger_z;
        scalar z_start = -0.5 * n_z * (module_ly - stagger_z);
        for (unsigned int iz = 0; iz < n_z; ++iz)
        {
            scalar z_pos = z_start + iz * z_step;
            for (unsigned int iphi = 0; iphi < n_phi; ++iphi)
            {
                scalar phi = -M_PI + iphi * phi_step;
                scalar r_addon = (iz % 2) ? -stagger_r : stagger_r;
                point3 p = {(r+r_addon) * std::cos(phi), (r+r_addon) * std::sin(phi), z_pos};
                vector3 z = {std::cos(phi + tilt_phi), std::sin(phi + tilt_phi), 0.};
                vector3 x = {z[1], -z[0], 0.};
                transforms.push_back(transform3(p, z, x));
            }
        }
        return {rectangle_bounds, transforms};
    }

} // namespace detray
