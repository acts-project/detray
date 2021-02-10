/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "core/surface.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/serializer2.hpp"
#include "grids/populator.hpp"
#include "utils/containers.hpp"
#include "utils/indexing.hpp"
#include "utils/enumerate.hpp"
#include "tools/local_object_finder.hpp"

#include <functional>

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
     * @returns a tuple for rectangle descriptions and transforms
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

    using barrel_point2 = __plugin::cylindrical2::point2;
    using endcap_point2 = __plugin::polar2::point2;
    using barrel_surface_finder = std::function<dvector<dindex>(const barrel_point2 &)>;

    /** This method creates a barrel description of surfaces
     * 
     * @param r The radius of the barrel cylinder
     * @param stagger_r The (optional) staggering in r
     * @param n_phi The number of modules in phi
     * @param tilt_phi The tilt of the modules in phi
     * @param overlap_rphi The overlap in rphi coverage (approximated)
     * @param barrel_z The length of the (sensitive) barrel in z
     * @param stagger_z The (optional) staggering in z
     * @param volume_inner_r The volume inner r
     * @param volume_outer_r The volume outer r
     * @param volume_half_z The volume half length
     * 
     * @returns a tuple for rectangle descriptions, transforms, object finders
     */
    dtuple<darray<scalar, 2>, dvector<transform3>, dvector<barrel_surface_finder>> 
            create_barrel_components(scalar r,
                                    scalar stagger_r,
                                    unsigned int n_phi,
                                    scalar tilt_phi,
                                    scalar overlap_rphi,
                                    scalar barrel_z,
                                    scalar stagger_z,
                                    unsigned int n_z,
                                    scalar volume_inner_r,
                                    scalar volume_outer_r,
                                    scalar volume_half_z)
    {
        // Estimate module dimensions
        scalar module_lx = 2 * r * M_PI * (1 + overlap_rphi) / n_phi;
        scalar module_ly = (barrel_z + (n_z - 1) * stagger_z) / n_z;
        darray<scalar, 2> rectangle_bounds = {0.5 * module_lx, 0.5 * module_ly};

        // Prepare the local finders 
        replace_populator<> replacer;
        serializer2 serializer;

        // The detector transforms
        dvector<transform3> transforms;
        scalar phi_step = 2 * M_PI / n_phi;
        scalar z_step = module_ly - stagger_z;
        scalar z_start = -0.5 * n_z * (module_ly - stagger_z);

        // Declare the inner, outer, ecn, ecp object finder 
        axis::circular<> rphi_axis_inner = {n_phi, -volume_inner_r*(M_PI+0.5*phi_step), volume_inner_r*(M_PI-0.5*phi_step) };
        axis::closed<> z_axis_inner = {n_z, -0.5 * barrel_z, 0.5 * barrel_z};
        axis::circular<> rphi_axis_outer = {n_phi, -volume_outer_r*(M_PI+0.5*phi_step), volume_outer_r*(M_PI-0.5*phi_step) };
        axis::closed<> z_axis_outer = {n_z, -0.5 * barrel_z, 0.5 * barrel_z};
        axis::closed<> r_axis_ecn = {1, volume_inner_r, volume_outer_r};
        axis::circular<> phi_axis_ecn = {n_phi, -M_PI-0.5*phi_step, M_PI-0.5*phi_step };
        axis::closed<> r_axis_ecp = {1, volume_inner_r, volume_outer_r};
        axis::circular<> phi_axis_ecp = {n_phi, -M_PI-0.5*phi_step, M_PI-0.5*phi_step };

        using barrel_grid = grid2<decltype(replacer), decltype(rphi_axis_inner), decltype(z_axis_inner), decltype(serializer)>;
        using ec_grid = grid2<decltype(replacer), decltype(r_axis_ecn), decltype(phi_axis_ecn), decltype(serializer)>;

        barrel_grid barrel_grid_inner(std::move(rphi_axis_inner), std::move(z_axis_inner));
        barrel_grid barrel_grid_outer(std::move(rphi_axis_outer), std::move(z_axis_inner));
        ec_grid ec_grid_n(std::move(r_axis_ecn), std::move(phi_axis_ecn));
        ec_grid ec_grid_p(std::move(r_axis_ecp), std::move(phi_axis_ecp));

        for (unsigned int iz = 0; iz < n_z; ++iz)
        {
            scalar z_pos = z_start + iz * z_step;
            for (unsigned int iphi = 0; iphi < n_phi; ++iphi)
            {
                
                scalar phi = -M_PI + iphi * phi_step;
                scalar cos_phi = std::cos(phi);

                barrel_grid_inner.populate(barrel_point2{volume_inner_r*cos_phi, z_pos}, transforms.size());
                barrel_grid_outer.populate(barrel_point2{volume_outer_r*cos_phi, z_pos}, transforms.size());
                if (iz == 0){
                    ec_grid_n.populate(endcap_point2{r, phi}, transforms.size());
                }
                if (iz == n_z-1){
                    ec_grid_p.populate(endcap_point2{r, phi}, transforms.size());
                }

                // The transform
                scalar r_addon = (iz % 2) ? -stagger_r : stagger_r;
                point3 p = {(r+r_addon) * std::cos(phi), (r+r_addon) * std::sin(phi), z_pos};
                vector3 z = {std::cos(phi + tilt_phi), std::sin(phi + tilt_phi), 0.};
                vector3 x = {z[1], -z[0], 0.};
                transforms.push_back(transform3(p, z, x));
            }
        }

        local_zone_finder<barrel_grid> inner_finder(std::move(barrel_grid_inner));
        local_zone_finder<barrel_grid> outer_finder(std::move(barrel_grid_outer));
        local_zone_finder<ec_grid> ecn_finder(std::move(ec_grid_n));
        local_zone_finder<ec_grid> ecp_finder(std::move(ec_grid_p));

        return {rectangle_bounds, transforms, {inner_finder, outer_finder, ecn_finder, ecp_finder} };
    }

} // namespace detray
