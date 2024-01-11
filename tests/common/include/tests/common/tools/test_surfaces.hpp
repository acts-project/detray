/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/grids/axis.hpp"
#include "detray/grids/grid2.hpp"
#include "detray/grids/populator.hpp"
#include "detray/grids/serializer2.hpp"
#include "detray/masks/rectangle2D.hpp"
#include "detray/tools/local_object_finder.hpp"
#include "detray/utils/ranges.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <functional>

namespace detray {

namespace {

vecmem::host_memory_resource host_mr;

// TODO: Remove Cyclic dependendy with benchmark_intersec_surfaces.inl types
enum plane_mask_ids : unsigned int {
    e_plane_rectangle2 = 0u,
};

enum plane_material_ids : unsigned int {
    e_plane_slab = 0u,
};

using transform3 = __plugin::transform3<detray::scalar>;
using point3 = __plugin::point3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;

using plane_mask_link_t = dtyped_index<plane_mask_ids, dindex>;
using plane_material_link_t = dtyped_index<plane_material_ids, dindex>;

using binned_neighborhood = darray<darray<dindex, 2>, 2>;

/// This method creates a number (distances.size()) planes along a direction
[[maybe_unused]] dvector<
    surface_descriptor<plane_mask_link_t, plane_material_link_t, transform3>>
planes_along_direction(const dvector<scalar> &distances, vector3 direction) {
    // Rotation matrix
    vector3 z = direction;
    vector3 x = vector::normalize(vector3{0.f, -z[2], z[1]});

    dvector<surface_descriptor<plane_mask_link_t, plane_material_link_t,
                               transform3>>
        surfaces;
    surfaces.reserve(distances.size());
    for (const auto [idx, d] : detray::views::enumerate(distances)) {
        vector3 t = d * direction;
        transform3 trf(t, z, x);
        plane_mask_link_t mask_link{plane_mask_ids::e_plane_rectangle2, idx};
        plane_material_link_t material_link{plane_material_ids::e_plane_slab,
                                            0u};
        surfaces.emplace_back(std::move(trf), std::move(mask_link),
                              std::move(material_link), 0u,
                              surface_id::e_sensitive);
        surfaces.back().set_index(idx);
    }
    return surfaces;
}

using cylinder_point2 = __plugin::point2<detray::scalar>;
using disc_point2 = __plugin::point2<detray::scalar>;
using endcap_surface_finder = std::function<dvector<dindex>(
    const disc_point2 &, const binned_neighborhood &)>;

/// This method creates an endcap description of surfaces
///
/// @param inner_r The inner radius of the endcap disc
/// @param outer_r The outer radius of the endcap disc
/// @param pos_z The nominal z position of the endcap disc
/// @param stagger_z The possible staggering in z
/// @param n_phi The number of modules in phi
/// @param overlap_phi The overlap in phi coverage (approximated)
/// @param volume_inner_r the volume inner radius for the local finder
/// @param volume_outer_r the volume outer radius for the local finder
/// @param volumer_min_z the volume minimum z
/// @param volume_max_z the volume maximumz
/// @param transform_offset The offset for the transform grid
///
/// @returns a tuple for rectangle descriptions and transforms
[[maybe_unused]] dtuple<darray<scalar, 3>, dvector<transform3>,
                        dvector<endcap_surface_finder>>
create_endcap_components(scalar inner_r, scalar outer_r, scalar pos_z,
                         scalar stagger_z, unsigned int n_phi,
                         scalar overlap_phi, scalar volume_inner_r,
                         scalar volume_outer_r, scalar volume_min_z,
                         scalar volume_max_z,
                         unsigned int transform_offset = 0u) {
    scalar module_inner_lx{2.f * inner_r * constant<scalar>::pi *
                           (1.f + overlap_phi) / static_cast<scalar>(n_phi)};
    scalar module_outer_lx{2.f * outer_r * constant<scalar>::pi *
                           (1.f + overlap_phi) / static_cast<scalar>(n_phi)};
    scalar module_hy{0.5f * (outer_r - inner_r)};

    darray<scalar, 3> trapezoid_values = {0.5f * module_inner_lx,
                                          0.5f * module_outer_lx, module_hy};
    scalar step_phi{2.f * constant<scalar>::pi / static_cast<scalar>(n_phi)};
    dvector<transform3> transforms;

    // Prepare the local finders
    serializer2 serializer;

    // Declare the inner, outer, ecn, ecp object finder

    using cylinder_grid = grid2<replace_populator, axis::circular,
                                axis::regular, decltype(serializer)>;
    using disc_grid = grid2<replace_populator, axis::regular, axis::circular,
                            decltype(serializer)>;

    typename cylinder_grid::axis_p0_type rphi_axis_inner = {
        n_phi, -volume_inner_r * (constant<scalar>::pi + 0.5f * step_phi),
        volume_inner_r * (constant<scalar>::pi - 0.5f * step_phi), host_mr};
    typename cylinder_grid::axis_p1_type z_axis_inner = {1u, volume_min_z,
                                                         volume_max_z, host_mr};
    typename cylinder_grid::axis_p0_type rphi_axis_outer = {
        n_phi, -volume_outer_r * (constant<scalar>::pi + 0.5f * step_phi),
        volume_outer_r * (constant<scalar>::pi - 0.5f * step_phi), host_mr};

    typename disc_grid::axis_p0_type r_axis_ecn = {1u, volume_inner_r,
                                                   volume_outer_r, host_mr};
    typename disc_grid::axis_p1_type phi_axis_ecn = {
        n_phi, -constant<scalar>::pi - 0.5f * step_phi,
        constant<scalar>::pi - 0.5f * step_phi, host_mr};
    typename disc_grid::axis_p0_type r_axis_ecp = {1u, volume_inner_r,
                                                   volume_outer_r, host_mr};
    typename disc_grid::axis_p1_type phi_axis_ecp = {
        n_phi, -constant<scalar>::pi - 0.5f * step_phi,
        constant<scalar>::pi - 0.5f * step_phi, host_mr};

    cylinder_grid ec_grid_inner(std::move(rphi_axis_inner),
                                std::move(z_axis_inner), host_mr);
    cylinder_grid ec_grid_outer(std::move(rphi_axis_outer),
                                std::move(z_axis_inner), host_mr);
    disc_grid ec_grid_n(std::move(r_axis_ecn), std::move(phi_axis_ecn),
                        host_mr);
    disc_grid ec_grid_p(std::move(r_axis_ecp), std::move(phi_axis_ecp),
                        host_mr);

    scalar r{0.5f * (inner_r + outer_r)};

    for (unsigned int iphi = 0u; iphi < n_phi; ++iphi) {
        scalar phi{-constant<scalar>::pi +
                   static_cast<scalar>(iphi) * step_phi};

        // Populate the grids
        dindex transform_idx{static_cast<dindex>(transforms.size()) +
                             transform_offset};
        ec_grid_inner.populate(cylinder_point2{volume_inner_r * phi, pos_z},
                               transform_idx + 0u);
        ec_grid_outer.populate(cylinder_point2{volume_outer_r * phi, pos_z},
                               transform_idx + 0u);
        ec_grid_n.populate(disc_point2{r, phi}, transform_idx + 0u);
        ec_grid_p.populate(disc_point2{r, phi}, transform_idx + 0u);

        scalar z_addon = (iphi % 2u) ? -stagger_z : stagger_z;
        scalar cos_phi = std::cos(phi);
        scalar sin_phi = std::sin(phi);
        point3 p = {r * cos_phi, r * sin_phi, pos_z + z_addon};
        vector3 z = {0.f, 0.f, 1.f};
        vector3 x = {sin_phi, -cos_phi, 0.f};
        transforms.push_back(transform3(p, z, x));
    }

    local_zone_finder<cylinder_grid> inner_finder(std::move(ec_grid_inner));
    local_zone_finder<cylinder_grid> outer_finder(std::move(ec_grid_outer));
    local_zone_finder<disc_grid> ecn_finder(std::move(ec_grid_n));
    local_zone_finder<disc_grid> ecp_finder(std::move(ec_grid_p));

    return {trapezoid_values, transforms,
            dvector<endcap_surface_finder>{inner_finder, outer_finder,
                                           ecn_finder, ecp_finder}};
}

using barrel_surface_finder = std::function<dvector<dindex>(
    const cylinder_point2 &, const binned_neighborhood &)>;

/// This method creates a barrel description of surfaces
///
/// @param r The radius of the barrel cylinder
/// @param stagger_r The (optional) staggering in r
/// @param n_phi The number of modules in phi
/// @param tilt_phi The tilt of the modules in phi
/// @param overlap_rphi The overlap in rphi coverage (approximated)
/// @param length_z The length of the (sensitive) barrel in z
/// @param overlap_z The (optional) staggering in z
/// @param volume_inner_r The volume inner r
/// @param volume_outer_r The volume outer r
/// @param volume_length_z The volume half length
/// @param transform_offset The offset for the transform grid
///
/// @returns a tuple for rectangle descriptions, transforms, object finders
[[maybe_unused]] dtuple<darray<scalar, 2>, dvector<transform3>,
                        dvector<barrel_surface_finder>>
create_barrel_components(scalar r, scalar stagger_r, unsigned int n_phi,
                         scalar tilt_phi, scalar overlap_rphi, scalar length_z,
                         scalar overlap_z, unsigned int n_z,
                         scalar volume_inner_r, scalar volume_outer_r,
                         scalar /*volume_half_z*/,
                         unsigned int transform_offset = 0u) {
    // Estimate module dimensions
    scalar module_lx{2.f * r * constant<scalar>::pi * (1.f + overlap_rphi) /
                     static_cast<scalar>(n_phi)};
    scalar module_ly{(length_z + static_cast<scalar>(n_z - 1u) * overlap_z) /
                     static_cast<scalar>(n_z)};
    darray<scalar, 2> rectangle_bounds = {0.5f * module_lx, 0.5f * module_ly};

    // Prepare the local finders
    serializer2 serializer;

    // The detector transforms
    dvector<transform3> transforms;
    scalar step_phi{2.f * constant<scalar>::pi / static_cast<scalar>(n_phi)};
    scalar step_z{module_ly - overlap_z};
    scalar start_z{-0.5f * static_cast<scalar>(n_z - 1u) *
                   (module_ly - overlap_z)};

    // Declare the inner, outer, ecn, ecp object finder
    using cylinder_grid = grid2<replace_populator, axis::circular,
                                axis::regular, decltype(serializer)>;
    using disc_grid = grid2<replace_populator, axis::regular, axis::circular,
                            decltype(serializer)>;

    typename cylinder_grid::axis_p0_type rphi_axis_inner = {
        n_phi, -volume_inner_r * (constant<scalar>::pi + 0.5f * step_phi),
        volume_inner_r * (constant<scalar>::pi - 0.5f * step_phi), host_mr};
    typename cylinder_grid::axis_p1_type z_axis_inner = {
        n_z, -0.5f * length_z, 0.5f * length_z, host_mr};
    typename cylinder_grid::axis_p0_type rphi_axis_outer = {
        n_phi, -volume_outer_r * (constant<scalar>::pi + 0.5f * step_phi),
        volume_outer_r * (constant<scalar>::pi - 0.5f * step_phi), host_mr};
    // axis::regular<> z_axis_outer = {n_z, -0.5f * length_z, 0.5f * length_z};
    typename disc_grid::axis_p0_type r_axis_ecn = {1, volume_inner_r,
                                                   volume_outer_r, host_mr};
    typename disc_grid::axis_p1_type phi_axis_ecn = {
        n_phi, -constant<scalar>::pi - 0.5f * step_phi,
        constant<scalar>::pi - 0.5f * step_phi, host_mr};
    typename disc_grid::axis_p0_type r_axis_ecp = {1, volume_inner_r,
                                                   volume_outer_r, host_mr};
    typename disc_grid::axis_p1_type phi_axis_ecp = {
        n_phi, -constant<scalar>::pi - 0.5f * step_phi,
        constant<scalar>::pi - 0.5f * step_phi, host_mr};

    cylinder_grid barrel_grid_inner(std::move(rphi_axis_inner),
                                    std::move(z_axis_inner), host_mr);
    cylinder_grid barrel_grid_outer(std::move(rphi_axis_outer),
                                    std::move(z_axis_inner), host_mr);
    disc_grid barrel_grid_n(std::move(r_axis_ecn), std::move(phi_axis_ecn),
                            host_mr);
    disc_grid barrel_grid_p(std::move(r_axis_ecp), std::move(phi_axis_ecp),
                            host_mr);

    for (unsigned int iz = 0u; iz < n_z; ++iz) {
        scalar pos_z = start_z + static_cast<scalar>(iz) * step_z;
        for (unsigned int iphi = 0u; iphi < n_phi; ++iphi) {
            scalar phi{-constant<scalar>::pi +
                       static_cast<scalar>(iphi) * step_phi};
            // Populate the grids
            dindex transform_idx{static_cast<dindex>(transforms.size()) +
                                 transform_offset};
            barrel_grid_inner.populate(
                cylinder_point2{volume_inner_r * phi, pos_z},
                transform_idx + 0u);
            barrel_grid_outer.populate(
                cylinder_point2{volume_outer_r * phi, pos_z},
                transform_idx + 0u);
            if (iz == 0u) {
                barrel_grid_n.populate(disc_point2{r, phi}, transform_idx + 0u);
            }
            if (iz == n_z - 1u) {
                barrel_grid_p.populate(disc_point2{r, phi}, transform_idx + 0u);
            }
            // Finally create the transform
            scalar r_addon = (iz % 2u) ? -stagger_r : stagger_r;
            point3 p = {(r + r_addon) * std::cos(phi),
                        (r + r_addon) * std::sin(phi), pos_z};
            vector3 z = {std::cos(phi + tilt_phi), std::sin(phi + tilt_phi),
                         0.f};
            vector3 x = {z[1], -z[0], 0.f};
            transforms.push_back(transform3(p, z, x));
        }
    }

    local_zone_finder<cylinder_grid> inner_finder(std::move(barrel_grid_inner));
    local_zone_finder<cylinder_grid> outer_finder(std::move(barrel_grid_outer));
    local_zone_finder<disc_grid> ecn_finder(std::move(barrel_grid_n));
    local_zone_finder<disc_grid> ecp_finder(std::move(barrel_grid_p));

    return {rectangle_bounds, transforms,
            dvector<barrel_surface_finder>{inner_finder, outer_finder,
                                           ecn_finder, ecp_finder}};
}

}  // anonymous namespace

}  // namespace detray
