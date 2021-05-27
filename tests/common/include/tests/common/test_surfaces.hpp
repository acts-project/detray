/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "core/surface_base.hpp"
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
    using point3 = point3;
    using vector3 = vector3;

    /** This method creates a number (distances.size()) planes along a direction 
    */
    dvector<surface_base<transform3>> planes_along_direction(dvector<scalar> distances, vector3 direction)
    {
        // Rotation matrix
        vector3 z = direction;
        vector3 x = normalize(vector3{0, -z[2], z[1]});

        dvector<surface_base<transform3>> return_surfaces;
        return_surfaces.reserve(distances.size());
        for (auto &d : distances)
        {
            vector3 t = d * direction;
            transform3 trf(t, z, x);
            return_surfaces.push_back(surface_base<transform3>{std::move(trf), 0, 0, false});
        }
        return return_surfaces;
    }

    using cylinder_point2 = __plugin::point2;
    using disc_point2 = __plugin::point2;
    using endcap_surface_finder = std::function<dvector<dindex>(const disc_point2 &, const darray<unsigned int, 2> &)>;

    /** This method creates a barrel description of surfaces
     * 
     * @param inner_r The inner radius of the endcap disc
     * @param outer_r The outer radius of the endcap disc
     * @param pos_z The nominal z position of the endcap disc
     * @param stagger_z The possible staggering in z 
     * @param n_phi The number of modules in phi
     * @param overlap_phi The overlap in phi coverage (approximated)
     * @param volume_inner_r the volume inner radius for the local finder
     * @param volume_outer_r the volume outer radius for the local finder
     * @param volumer_min_z the volume minimum z
     * @param volume_max_z the volume maximumz
     * @param transform_offset The offset for the transform grid
     * 
     * @returns a tuple for rectangle descriptions and transforms
     */
    dtuple<darray<scalar, 3>, dvector<transform3>, dvector<endcap_surface_finder>>
    create_endcap_components(scalar inner_r,
                             scalar outer_r,
                             scalar pos_z,
                             scalar stagger_z,
                             unsigned int n_phi,
                             scalar overlap_phi,
                             scalar volume_inner_r,
                             scalar volume_outer_r,
                             scalar volume_min_z,
                             scalar volume_max_z,
                             unsigned int transform_offset = 0)
    {
        scalar module_inner_lx = 2. * inner_r * M_PI * (1. + overlap_phi) / (n_phi);
        scalar module_outer_lx = 2. * outer_r * M_PI * (1. + overlap_phi) / (n_phi);
        scalar module_hy = 0.5 * (outer_r - inner_r);

        darray<scalar, 3> trapezoid_values = {static_cast<scalar>(0.5 * module_inner_lx), static_cast<scalar>(0.5 * module_outer_lx), module_hy};
        scalar step_phi = 2 * M_PI / n_phi;
        dvector<transform3> transforms;

        // Prepare the local finders
        replace_populator<> replacer;
        serializer2 serializer;

        // Declare the inner, outer, ecn, ecp object finder
        axis::circular rphi_axis_inner = {n_phi,
                                          static_cast<scalar>(-volume_inner_r * (M_PI + 0.5 * step_phi)),
                                          static_cast<scalar>(volume_inner_r * (M_PI - 0.5 * step_phi))};
        axis::regular z_axis_inner = {1, volume_min_z, volume_max_z};
        axis::circular rphi_axis_outer = {n_phi,
                                          static_cast<scalar>(-volume_outer_r * (M_PI + 0.5 * step_phi)),
                                          static_cast<scalar>(volume_outer_r * (M_PI - 0.5 * step_phi))};
        axis::regular z_axis_outer = {1, volume_min_z, volume_max_z};
        axis::regular r_axis_ecn = {1, volume_inner_r, volume_outer_r};
        axis::circular phi_axis_ecn = {n_phi,
                                       static_cast<scalar>(-M_PI - 0.5 * step_phi),
                                       static_cast<scalar>(M_PI - 0.5 * step_phi)};
        axis::regular r_axis_ecp = {1, volume_inner_r, volume_outer_r};
        axis::circular phi_axis_ecp = {n_phi,
                                       static_cast<scalar>(-M_PI - 0.5 * step_phi),
                                       static_cast<scalar>(M_PI - 0.5 * step_phi)};

        using cylinder_grid = grid2<decltype(replacer), decltype(rphi_axis_inner), decltype(z_axis_inner), decltype(serializer)>;
        using disc_grid = grid2<decltype(replacer), decltype(r_axis_ecn), decltype(phi_axis_ecn), decltype(serializer)>;

        cylinder_grid ec_grid_inner(std::move(rphi_axis_inner), std::move(z_axis_inner));
        cylinder_grid ec_grid_outer(std::move(rphi_axis_outer), std::move(z_axis_inner));
        disc_grid ec_grid_n(std::move(r_axis_ecn), std::move(phi_axis_ecn));
        disc_grid ec_grid_p(std::move(r_axis_ecp), std::move(phi_axis_ecp));

        scalar r = 0.5 * (inner_r + outer_r);

        for (unsigned int iphi = 0; iphi < n_phi; ++iphi)
        {
            scalar phi = -M_PI + iphi * step_phi;

            // Populate the grids
            ec_grid_inner.populate(cylinder_point2{volume_inner_r * phi, pos_z}, transforms.size() + transform_offset);
            ec_grid_outer.populate(cylinder_point2{volume_outer_r * phi, pos_z}, transforms.size() + transform_offset);
            ec_grid_n.populate(disc_point2{r, phi}, transforms.size() + transform_offset);
            ec_grid_p.populate(disc_point2{r, phi}, transforms.size() + transform_offset);

            scalar z_addon = (iphi % 2) ? -stagger_z : stagger_z;
            scalar cos_phi = std::cos(phi);
            scalar sin_phi = std::sin(phi);
            point3 p = {r * cos_phi, r * sin_phi, pos_z + z_addon};
            vector3 z = {0., 0., 1.};
            vector3 x = {sin_phi, -cos_phi, 0.};
            transforms.push_back(transform3(p, z, x));
        }

        local_zone_finder<cylinder_grid> inner_finder(std::move(ec_grid_inner));
        local_zone_finder<cylinder_grid> outer_finder(std::move(ec_grid_outer));
        local_zone_finder<disc_grid> ecn_finder(std::move(ec_grid_n));
        local_zone_finder<disc_grid> ecp_finder(std::move(ec_grid_p));

        return {trapezoid_values, transforms, {inner_finder, outer_finder, ecn_finder, ecp_finder}};
    }

    using barrel_surface_finder = std::function<dvector<dindex>(const cylinder_point2 &, const darray<unsigned int, 2> &)>;

    /** This method creates a barrel description of surfaces
     * 
     * @param r The radius of the barrel cylinder
     * @param stagger_r The (optional) staggering in r
     * @param n_phi The number of modules in phi
     * @param tilt_phi The tilt of the modules in phi
     * @param overlap_rphi The overlap in rphi coverage (approximated)
     * @param length_z The length of the (sensitive) barrel in z
     * @param overlap_z The (optional) staggering in z
     * @param volume_inner_r The volume inner r
     * @param volume_outer_r The volume outer r
     * @param volume_length_z The volume half length
     * @param transform_offset The offset for the transform grid
     * 
     * @returns a tuple for rectangle descriptions, transforms, object finders
     */
    dtuple<darray<scalar, 2>, dvector<transform3>, dvector<barrel_surface_finder>>
    create_barrel_components(scalar r,
                             scalar stagger_r,
                             unsigned int n_phi,
                             scalar tilt_phi,
                             scalar overlap_rphi,
                             scalar length_z,
                             scalar overlap_z,
                             unsigned int n_z,
                             scalar volume_inner_r,
                             scalar volume_outer_r,
                             scalar volume_half_z,
                             unsigned int transform_offset = 0)
    {
        // Estimate module dimensions
        scalar module_lx = 2 * r * M_PI * (1 + overlap_rphi) / n_phi;
        scalar module_ly = (length_z + (n_z - 1) * overlap_z) / n_z;
        darray<scalar, 2> rectangle_bounds = {static_cast<scalar>(0.5 * module_lx), static_cast<scalar>(0.5 * module_ly)};

        // Prepare the local finders
        replace_populator<> replacer;
        serializer2 serializer;

        // The detector transforms
        dvector<transform3> transforms;
        scalar step_phi = 2 * M_PI / n_phi;
        scalar step_z = module_ly - overlap_z;
        scalar start_z = -0.5 * (n_z - 1) * (module_ly - overlap_z);

        // Declare the inner, outer, ecn, ecp object finder
        axis::circular rphi_axis_inner = {n_phi,
                                          static_cast<scalar>(-volume_inner_r * (M_PI + 0.5 * step_phi)),
                                          static_cast<scalar>(volume_inner_r * (M_PI - 0.5 * step_phi))};
        axis::regular z_axis_inner = {n_z, static_cast<scalar>(-0.5 * length_z), static_cast<scalar>(0.5 * length_z)};
        axis::circular rphi_axis_outer = {n_phi,
                                          static_cast<scalar>(-volume_outer_r * (M_PI + 0.5 * step_phi)),
                                          static_cast<scalar>(volume_outer_r * (M_PI - 0.5 * step_phi))};
        axis::regular z_axis_outer = {n_z, static_cast<scalar>(-0.5 * length_z), static_cast<scalar>(0.5 * length_z)};
        axis::regular r_axis_ecn = {1, volume_inner_r, volume_outer_r};
        axis::circular phi_axis_ecn = {n_phi,
                                       static_cast<scalar>(-M_PI - 0.5 * step_phi),
                                       static_cast<scalar>(M_PI - 0.5 * step_phi)};
        axis::regular r_axis_ecp = {1, volume_inner_r, volume_outer_r};
        axis::circular phi_axis_ecp = {n_phi,
                                       static_cast<scalar>(-M_PI - 0.5 * step_phi),
                                       static_cast<scalar>(M_PI - 0.5 * step_phi)};

        using cylinder_grid = grid2<decltype(replacer), decltype(rphi_axis_inner), decltype(z_axis_inner), decltype(serializer)>;
        using disc_grid = grid2<decltype(replacer), decltype(r_axis_ecn), decltype(phi_axis_ecn), decltype(serializer)>;

        cylinder_grid barrel_grid_inner(std::move(rphi_axis_inner), std::move(z_axis_inner));
        cylinder_grid barrel_grid_outer(std::move(rphi_axis_outer), std::move(z_axis_inner));
        disc_grid barrel_grid_n(std::move(r_axis_ecn), std::move(phi_axis_ecn));
        disc_grid barrel_grid_p(std::move(r_axis_ecp), std::move(phi_axis_ecp));

        for (unsigned int iz = 0; iz < n_z; ++iz)
        {
            scalar pos_z = start_z + iz * step_z;
            for (unsigned int iphi = 0; iphi < n_phi; ++iphi)
            {
                scalar phi = -M_PI + iphi * step_phi;
                // Populate the grids
                barrel_grid_inner.populate(cylinder_point2{volume_inner_r * phi, pos_z}, transforms.size() + transform_offset);
                barrel_grid_outer.populate(cylinder_point2{volume_outer_r * phi, pos_z}, transforms.size() + transform_offset);
                if (iz == 0)
                {
                    barrel_grid_n.populate(disc_point2{r, phi}, transforms.size() + transform_offset);
                }
                if (iz == n_z - 1)
                {
                    barrel_grid_p.populate(disc_point2{r, phi}, transforms.size() + transform_offset);
                }
                // Finally create the transform
                scalar r_addon = (iz % 2) ? -stagger_r : stagger_r;
                point3 p = {(r + r_addon) * std::cos(phi), (r + r_addon) * std::sin(phi), pos_z};
                vector3 z = {std::cos(phi + tilt_phi), std::sin(phi + tilt_phi), 0.};
                vector3 x = {z[1], -z[0], 0.};
                transforms.push_back(transform3(p, z, x));
            }
        }

        local_zone_finder<cylinder_grid> inner_finder(std::move(barrel_grid_inner));
        local_zone_finder<cylinder_grid> outer_finder(std::move(barrel_grid_outer));
        local_zone_finder<disc_grid> ecn_finder(std::move(barrel_grid_n));
        local_zone_finder<disc_grid> ecp_finder(std::move(barrel_grid_p));

        return {rectangle_bounds, transforms, {inner_finder, outer_finder, ecn_finder, ecp_finder}};
    }

} // namespace detray
