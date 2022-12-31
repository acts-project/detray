/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/detector_metadata.hpp"
#include "detray/materials/predefined_materials.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/field.hpp>

// System include(s)
#include <climits>
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace detray {

namespace {

using point3 = __plugin::point3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;
using point2 = __plugin::point2<detray::scalar>;

/** Function that adds a cylinder portal.
 *
 * @tparam cylinder_id default cylinder id
 *
 * @param volume_idx volume the portal should be added to
 * @param ctx geometric context
 * @param surfaces container to add new surface to
 * @param masks container to add new cylinder mask to
 * @param transforms container to add new transform to
 * @param r raius of the cylinder
 * @param lower_z lower extend in z
 * @param upper_z upper extend in z
 * @param volume_link link to next volume for the masks
 */
template <typename context_t, typename surface_container_t,
          typename mask_container_t, typename material_container_t,
          typename transform_container_t, typename volume_links>
inline void add_cylinder_surface(
    const dindex volume_idx, context_t &ctx, surface_container_t &surfaces,
    mask_container_t &masks, material_container_t &materials,
    transform_container_t &transforms, const scalar r, const scalar lower_z,
    const scalar upper_z, const volume_links volume_link,
    const material<scalar> &mat, const scalar thickness) {
    using surface_type = typename surface_container_t::value_type;
    using mask_id = typename surface_type::mask_id;
    using mask_link_type = typename surface_type::mask_link;
    using material_id = typename surface_type::material_id;
    using material_link_type = typename surface_type::material_link;

    constexpr auto cylinder_id = mask_id::e_portal_cylinder2;
    constexpr auto slab_id = material_id::e_slab;

    const scalar min_z{std::min(lower_z, upper_z)};
    const scalar max_z{std::max(lower_z, upper_z)};

    // translation
    point3 tsl{0., 0., 0};

    // add transform and masks
    transforms.emplace_back(ctx, tsl);
    masks.template emplace_back<cylinder_id>(empty_context{}, volume_link, r,
                                             min_z, max_z);

    // Add material slab
    materials.template emplace_back<slab_id>(empty_context{}, mat, thickness);

    // add surface
    mask_link_type mask_link{cylinder_id,
                             masks.template size<cylinder_id>() - 1};
    material_link_type material_link{slab_id,
                                     materials.template size<slab_id>() - 1};
    const surface_id sf_id = (volume_link != volume_idx)
                                 ? surface_id::e_portal
                                 : surface_id::e_passive;
    surfaces.emplace_back(transforms.size(ctx) - 1, mask_link, material_link,
                          volume_idx, dindex_invalid, sf_id);
}

/** Function that adds a disc portal.
 *
 * @tparam disc_id default disc id
 *
 * @param volume_idx volume the portal should be added to
 * @param ctx geometric context
 * @param surfaces container to add new surface to
 * @param masks container to add new cylinder mask to
 * @param transforms container to add new transform to
 * @param min_r lower radius of the disc
 * @param max_r upper radius of the disc
 * @param z z position of the disc
 * @param volume_link link to next volume for the masks
 */
template <typename context_t, typename surface_container_t,
          typename mask_container_t, typename material_container_t,
          typename transform_container_t, typename volume_links>
inline void add_disc_surface(
    const dindex volume_idx, context_t &ctx, surface_container_t &surfaces,
    mask_container_t &masks, material_container_t &materials,
    transform_container_t &transforms, const scalar inner_r,
    const scalar outer_r, const scalar z, const volume_links volume_link,
    const material<scalar> &mat, const scalar thickness) {
    using surface_type = typename surface_container_t::value_type;
    using mask_id = typename surface_type::mask_id;
    using mask_link_type = typename surface_type::mask_link;
    using material_id = typename surface_type::material_id;
    using material_link_type = typename surface_type::material_link;

    constexpr auto disc_id = mask_id::e_portal_ring2;
    constexpr auto slab_id = material_id::e_slab;

    const scalar min_r{std::min(inner_r, outer_r)};
    const scalar max_r{std::max(inner_r, outer_r)};

    // translation
    point3 tsl{0., 0., z};

    // add transform and mask
    transforms.emplace_back(ctx, tsl);
    masks.template emplace_back<disc_id>(empty_context{}, volume_link, min_r,
                                         max_r);

    // Add material slab
    materials.template emplace_back<slab_id>(empty_context{}, mat, thickness);

    // add surface
    mask_link_type mask_link{disc_id, masks.template size<disc_id>() - 1};
    material_link_type material_link{slab_id,
                                     materials.template size<slab_id>() - 1};
    const surface_id sf_id = (volume_link != volume_idx)
                                 ? surface_id::e_portal
                                 : surface_id::e_sensitive;
    surfaces.emplace_back(transforms.size(ctx) - 1, mask_link, material_link,
                          volume_idx, dindex_invalid, sf_id);
}

/** Function that adds a generic cylinder volume, using a factory for contained
 *  module surfaces.
 *
 * @tparam detector_t the detector type
 * @tparam factory_t type of module factory. Must be callable on containers.
 *
 * @param det detector the volume should be added to
 * @param resource vecmem memory resource for the temporary containers
 * @param ctx geometry context
 * @param lay_inner_r inner radius of volume
 * @param lay_outer_r outer radius of volume
 * @param lay_neg_r lower extend of volume
 * @param lay_pos_r upper extend of volume
 * @param volume_links volume links for the portals of the volume
 * @param module_factory functor that adds module surfaces to volume
 */

template <
    typename detector_t, typename factory_t,
    std::enable_if_t<
        std::is_invocable_v<factory_t, typename detector_t::geometry_context &,
                            typename detector_t::volume_type &,
                            typename detector_t::surface_container_t &,
                            typename detector_t::mask_container &,
                            typename detector_t::material_container &,
                            typename detector_t::transform_container &>,
        bool> = true>
void create_cyl_volume(
    detector_t &det, vecmem::memory_resource &resource,
    typename detector_t::geometry_context &ctx, const scalar lay_inner_r,
    const scalar lay_outer_r, const scalar lay_neg_z, const scalar lay_pos_z,
    const std::vector<typename detector_t::surface_type::volume_link_type>
        &volume_links,
    factory_t &module_factory) {
    // volume bounds
    const scalar inner_r = std::min(lay_inner_r, lay_outer_r);
    const scalar outer_r = std::max(lay_inner_r, lay_outer_r);
    const scalar lower_z = std::min(lay_neg_z, lay_pos_z);
    const scalar upper_z = std::max(lay_neg_z, lay_pos_z);

    auto &cyl_volume =
        det.new_volume(volume_id::e_cylinder,
                       {inner_r, outer_r, lower_z, upper_z, -M_PI, M_PI});
    cyl_volume.set_link(detector_t::sf_finders::id::e_default, 0);

    // Add module surfaces to volume
    typename detector_t::surface_container_t surfaces(&resource);
    typename detector_t::mask_container masks(resource);
    typename detector_t::material_container materials(resource);
    typename detector_t::transform_container transforms(resource);

    // fill the surfaces
    module_factory(ctx, cyl_volume, surfaces, masks, materials, transforms);

    // negative and positive, inner and outer portal surface
    add_cylinder_surface(cyl_volume.index(), ctx, surfaces, masks, materials,
                         transforms, inner_r, lower_z, upper_z, volume_links[0],
                         vacuum<scalar>(), 0. * unit<scalar>::mm);
    add_cylinder_surface(cyl_volume.index(), ctx, surfaces, masks, materials,
                         transforms, outer_r, lower_z, upper_z, volume_links[1],
                         vacuum<scalar>(), 0. * unit<scalar>::mm);
    add_disc_surface(cyl_volume.index(), ctx, surfaces, masks, materials,
                     transforms, inner_r, outer_r, lower_z, volume_links[2],
                     vacuum<scalar>(), 0. * unit<scalar>::mm);
    add_disc_surface(cyl_volume.index(), ctx, surfaces, masks, materials,
                     transforms, inner_r, outer_r, upper_z, volume_links[3],
                     vacuum<scalar>(), 0. * unit<scalar>::mm);

    det.add_objects_per_volume(ctx, cyl_volume, surfaces, masks, transforms,
                               materials);
}

/** Helper function that creates a layer of rectangular barrel modules.
 *
 * @tparam rctangle_id default rectangle id
 *
 * @param ctx geometric context
 * @param volume_idx volume the portal should be added to
 * @param surfaces container to add new surface to
 * @param masks container to add new cylinder mask to
 * @param transforms container to add new transform to
 * @param cfg config struct for module creation
 */
template <typename context_t, typename volume_type,
          typename surface_container_t, typename mask_container_t,
          typename material_container_t, typename transform_container_t,
          typename config_t>
inline void create_barrel_modules(context_t &ctx, volume_type &vol,
                                  surface_container_t &surfaces,
                                  mask_container_t &masks,
                                  material_container_t &materials,
                                  transform_container_t &transforms,
                                  config_t cfg) {
    using surface_type = typename surface_container_t::value_type;
    using volume_link_t = typename surface_type::volume_link_type;
    using mask_id = typename surface_type::mask_id;
    using mask_link_type = typename surface_type::mask_link;
    using material_id = typename surface_type::material_id;
    using material_link_type = typename surface_type::material_link;

    constexpr auto rectangle_id = mask_id::e_rectangle2;
    constexpr auto slab_id = material_id::e_slab;

    auto volume_idx = vol.index();
    volume_link_t mask_volume_link{volume_idx};

    // Create the module centers

    // surface grid bins
    int n_phi_bins = cfg.m_binning.first;
    int n_z_bins = cfg.m_binning.second;
    // module positions
    detray::dvector<point3> m_centers;
    m_centers.reserve(n_phi_bins * n_z_bins);

    // prep work
    scalar pi{static_cast<scalar>(M_PI)};
    scalar phi_step{scalar{2} * pi / (n_phi_bins)};
    scalar min_phi{-pi + scalar{0.5} * phi_step};

    scalar z_start{scalar{-0.5} * (n_z_bins - 1) *
                   (scalar{2} * cfg.m_half_y - cfg.m_long_overlap)};
    scalar z_step{(std::abs(z_start) - z_start) / (n_z_bins - 1)};

    // loop over the z bins
    for (size_t z_bin = 0; z_bin < size_t(n_z_bins); ++z_bin) {
        // prepare z and r
        scalar m_z{z_start + z_bin * z_step};
        scalar m_r{(z_bin % 2) != 0u
                       ? cfg.layer_r - scalar{0.5} * cfg.m_radial_stagger
                       : cfg.layer_r + scalar{0.5} * cfg.m_radial_stagger};
        for (size_t phiBin = 0; phiBin < size_t(n_phi_bins); ++phiBin) {
            // calculate the current phi value
            scalar m_phi{min_phi + phiBin * phi_step};
            m_centers.push_back(
                point3{m_r * std::cos(m_phi), m_r * std::sin(m_phi), m_z});
        }
    }

    // Create geometry data
    for (auto &m_center : m_centers) {

        // Surfaces with the linking into the local containers
        mask_link_type mask_link = {rectangle_id,
                                    masks.template size<rectangle_id>()};
        material_link_type material_link{slab_id,
                                         materials.template size<slab_id>()};
        const auto trf_index = transforms.size(ctx);
        surfaces.emplace_back(trf_index, mask_link, material_link, volume_idx,
                              dindex_invalid, surface_id::e_sensitive);

        // The rectangle bounds for this module
        masks.template emplace_back<rectangle_id>(
            empty_context{}, mask_volume_link, cfg.m_half_x, cfg.m_half_y);
        materials.template emplace_back<slab_id>(empty_context{}, cfg.mat,
                                                 cfg.thickness);

        // Build the transform
        // The local phi
        scalar m_phi{algebra::getter::phi(m_center)};
        // Local z axis is the normal vector
        vector3 m_local_z{std::cos(m_phi + cfg.m_tilt_phi),
                          std::sin(m_phi + cfg.m_tilt_phi), 0.};
        // Local x axis the normal to local y,z
        vector3 m_local_x{-std::sin(m_phi + cfg.m_tilt_phi),
                          std::cos(m_phi + cfg.m_tilt_phi), 0.};

        // Create the module transform
        transforms.emplace_back(ctx, m_center, m_local_z, m_local_x);
    }
}

/** Helper function that creates a surface grid of rectangular barrel modules.
 *
 * @param surfaces_grid the grid to be created with proper axes
 * @param resource vecmem memory resource
 * @param cfg config struct for module creation
 */
template <typename detector_t, typename config_t>
inline void add_cylinder_grid(const typename detector_t::geometry_context &ctx,
                              typename detector_t::volume_type &vol,
                              detector_t &det, const config_t &cfg) {
    // Get relevant ids
    // using geo_obj_ids = typename detector_t::geo_obj_ids;

    constexpr auto cyl_id = detector_t::masks::id::e_portal_cylinder2;
    constexpr auto grid_id = detector_t::sf_finders::id::e_cylinder_grid;

    using cyl_grid_t =
        typename detector_t::surface_container::template get_type<grid_id>;
    auto gbuilder = grid_builder<detector_t, cyl_grid_t, detail::fill_by_pos>{};

    // The cylinder portals are at the end of the surface range by construction
    auto portal_mask_idx = (det.surfaces(vol).end() - 3)->mask().index();
    const auto &cyl_mask =
        det.mask_store().template get<cyl_id>().at(portal_mask_idx);

    // approximate binning for the barrel sensors
    std::size_t n_phi_bins{cfg.m_binning.first};
    std::size_t n_z_bins{cfg.m_binning.second};

    // Add new grid to the detector
    gbuilder.init_grid(cyl_mask, {n_phi_bins, n_z_bins});
    gbuilder.fill_grid(det, vol, ctx);

    det.surface_store().template push_back<grid_id>(gbuilder());
    // vol.set_link(grid_id,
    //                   det.surface_store().template size<grid_id>() - 1);
}

/** Helper function that creates a surface grid of trapezoidal endcap modules.
 *
 * @param vol the detector volume that should be equipped with a grid
 * @param surfaces_grid the grid to be created with proper axes
 * @param resource vecmem memory resource
 * @param cfg config struct for module creation
 */
template <typename detector_t, typename config_t>
inline void add_disc_grid(const typename detector_t::geometry_context &ctx,
                          typename detector_t::volume_type &vol,
                          detector_t &det, const config_t &cfg) {
    // Get relevant ids
    // using geo_obj_ids = typename detector_t::geo_obj_ids;

    constexpr auto disc_id = detector_t::masks::id::e_portal_ring2;
    constexpr auto grid_id = detector_t::sf_finders::id::e_disc_grid;

    using disc_grid_t =
        typename detector_t::surface_container::template get_type<grid_id>;
    auto gbuilder =
        grid_builder<detector_t, disc_grid_t, detray::detail::fill_by_pos>{};

    // The cylinder portals are at the end of the surface range by construction
    auto portal_mask_idx = det.surfaces(vol).back().mask().index();
    const auto &disc_mask =
        det.mask_store().template get<disc_id>().at(portal_mask_idx);

    // Add new grid to the detector
    gbuilder.init_grid(disc_mask,
                       {cfg.disc_binning.size(), cfg.disc_binning.front()});
    gbuilder.fill_grid(det, vol, ctx);

    det.surface_store().template push_back<grid_id>(gbuilder());
    // vol.set_link(grid_id,
    //                   det.surface_store().template size<grid_id>() - 1);
}

/** Helper method for positioning of modules in an endcap ring
 *
 * @param z is the z position of the ring
 * @param radius is the ring radius
 * @param phi_stagger is the radial staggering along phi
 * @param phi_sub_stagger is the overlap of the modules
 * @param n_phi_bins is the number of bins in phi
 *
 * @return a vector of the module positions in a ring
 */
inline auto module_positions_ring(scalar z, scalar radius, scalar phi_stagger,
                                  scalar phi_sub_stagger, int n_phi_bins) {
    // create and fill the positions
    std::vector<vector3> r_positions;
    r_positions.reserve(n_phi_bins);

    // prep work
    scalar pi{static_cast<scalar>(M_PI)};
    scalar phi_step{scalar{2} * pi / (n_phi_bins)};
    scalar min_phi{-pi + scalar{0.5} * phi_step};

    for (size_t iphi = 0; iphi < size_t(n_phi_bins); ++iphi) {
        // if we have a phi sub stagger presents
        scalar rzs = 0.;
        // phi stagger affects 0 vs 1, 2 vs 3 ... etc
        // -> only works if it is a %4
        // phi sub stagger affects 2 vs 4, 1 vs 3 etc.
        if (phi_sub_stagger != 0. && !(n_phi_bins % 4)) {
            // switch sides
            if (!(iphi % 4)) {
                rzs = phi_sub_stagger;
            } else if (!((iphi + 1) % 4)) {
                rzs = -phi_sub_stagger;
            }
        }
        // the module phi
        scalar phi{min_phi + iphi * phi_step};
        // main z position depending on phi bin
        scalar rz{iphi % 2 ? z - scalar{0.5} * phi_stagger
                           : z + scalar{0.5} * phi_stagger};
        r_positions.push_back(
            vector3{radius * std::cos(phi), radius * std::sin(phi), rz + rzs});
    }
    return r_positions;
}

/** Helper function that creates a layer of trapezoidal endcap modules.
 *
 * @tparam trapezoid_id default trapezoid id
 *
 * @param ctx geometric context
 * @param volume_idx volume the portal should be added to
 * @param surfaces container to add new surface to
 * @param masks container to add new cylinder mask to
 * @param transforms container to add new transform to
 * @param cfg config struct for module creation
 */
template <typename context_t, typename volume_type,
          typename surface_container_t, typename mask_container_t,
          typename material_container_t, typename transform_container_t,
          typename config_t>
void create_endcap_modules(context_t &ctx, volume_type &vol,
                           surface_container_t &surfaces,
                           mask_container_t &masks,
                           material_container_t &materials,
                           transform_container_t &transforms, config_t cfg) {
    using surface_type = typename surface_container_t::value_type;
    using volume_link_t = typename surface_type::volume_link_type;
    using mask_id = typename surface_type::mask_id;
    using mask_link_type = typename surface_type::mask_link;
    using material_id = typename surface_type::material_id;
    using material_link_type = typename surface_type::material_link;

    constexpr auto trapezoid_id = mask_id::e_trapezoid2;
    constexpr auto slab_id = material_id::e_slab;

    auto volume_idx = vol.index();
    volume_link_t mask_volume_link{volume_idx};

    // calculate the radii of the rings
    std::vector<scalar> radii;
    // calculate the radial borders
    // std::vector<scalar> radial_boarders;
    // the radial span of the disc
    scalar delta_r{cfg.outer_r - cfg.inner_r};

    // Only one ring
    if (cfg.disc_binning.size() == 1) {
        radii.push_back(scalar{0.5} * (cfg.inner_r + cfg.outer_r));
        // radial_boarders = {inner_r, outer_r};
    } else {
        // sum up the total length of the modules along r
        scalar tot_length{0};
        for (auto &m_hlength : cfg.m_half_y) {
            tot_length += 2 * m_hlength + 0.5;
        }
        // now calculate the overlap (equal pay)
        scalar r_overlap{(tot_length - delta_r) / (cfg.m_half_y.size() - 1)};
        // and now fill the radii and gaps
        scalar prev_r{cfg.inner_r};
        scalar prev_hl{0};
        scalar prev_ol{0};
        // remember the radial boarders
        // radial_boarders.push_back(inner_r);
        for (auto &m_hlength : cfg.m_half_y) {
            // calculate the radius
            radii.push_back(prev_r + prev_hl - prev_ol + m_hlength);
            prev_r = radii.back();
            prev_ol = r_overlap;
            prev_hl = m_hlength;
            // and register the radial boarder
            // radial_boarders.push_back(prev_r + scalar{2} * prev_hl -
            // scalar{0.5} * prev_ol);
        }
    }

    // now build the modules in every ring
    for (size_t ir = 0; ir < radii.size(); ++ir) {
        // generate the z value
        // convention inner ring is closer to origin : makes sense
        scalar rz{
            radii.size() == 1
                ? cfg.edc_position
                : (ir % 2 ? cfg.edc_position + scalar{0.5} * cfg.ring_stagger
                          : cfg.edc_position - scalar{0.5} * cfg.ring_stagger)};
        // fill the ring module positions
        scalar ps_stagger{
            cfg.m_phi_sub_stagger.size() ? cfg.m_phi_sub_stagger[ir] : 0};

        std::vector<point3> r_postitions =
            module_positions_ring(rz, radii[ir], cfg.m_phi_stagger[ir],
                                  ps_stagger, cfg.disc_binning[ir]);

        // Build the geometrical objects
        for (const auto &m_position : r_postitions) {
            // trapezoid mask
            mask_link_type mask_link{trapezoid_id,
                                     masks.template size<trapezoid_id>()};
            material_link_type material_link{
                slab_id, materials.template size<slab_id>()};

            masks.template emplace_back<trapezoid_id>(
                empty_context{}, mask_volume_link, cfg.m_half_x_min_y[ir],
                cfg.m_half_x_max_y[ir], cfg.m_half_y[ir],
                static_cast<scalar>(1. / (2. * cfg.m_half_y[ir])));

            materials.template emplace_back<slab_id>(empty_context{}, cfg.mat,
                                                     cfg.thickness);

            // Surfaces with the linking into the local containers
            surfaces.emplace_back(transforms.size(ctx), mask_link,
                                  material_link, volume_idx, dindex_invalid,
                                  surface_id::e_sensitive);

            // the module transform from the position
            scalar m_phi{algebra::getter::phi(m_position)};
            // the center position of the modules
            point3 m_center{static_cast<scalar>(cfg.side) * m_position};
            // the rotation matrix of the module
            vector3 m_local_y{std::cos(m_phi), std::sin(m_phi), 0.};
            // take different axis to have the same readout direction
            vector3 m_local_z{0., 0., cfg.side * scalar{1.}};
            vector3 m_local_x = algebra::vector::cross(m_local_y, m_local_z);

            // Create the module transform
            transforms.emplace_back(ctx, m_center, m_local_z, m_local_x);
        }
    }
}

/** Helper method for creating a beampipe with enclosing volume.
 *
 * @param det detector the volume should be added to
 * @param resource vecmem memory resource for the temporary containers
 * @param ctx geometry context
 * @param n_layers number of endcap layers that contain modules
 * @param edc_lay_sizes extend of the outer cylinder portal surfaces in z
 * @param beampipe_vol_size inner and outer radious of the beampipe volume
 * @param beampipe_r radius of the beampipe surface
 * @param brl_half_z half length of the barrel region in z
 */
template <typename detector_t>
inline void add_beampipe(
    detector_t &det, vecmem::memory_resource &resource,
    typename detector_t::geometry_context &ctx, const std::size_t n_edc_layers,
    const std::size_t n_brl_layers,
    const std::vector<std::pair<scalar, scalar>> &edc_lay_sizes,
    const std::pair<scalar, scalar> &beampipe_vol_size, const scalar beampipe_r,
    const scalar brl_half_z, const scalar edc_inner_r) {

    const dindex leaving_world = dindex_invalid;

    scalar max_z =
        n_edc_layers <= 0 ? brl_half_z : edc_lay_sizes[n_edc_layers - 1].second;
    scalar min_z = -max_z;

    typename detector_t::surface_container_t surfaces(&resource);
    typename detector_t::mask_container masks(resource);
    typename detector_t::material_container materials(resource);
    typename detector_t::transform_container transforms(resource);

    auto &beampipe =
        det.new_volume(volume_id::e_cylinder,
                       {beampipe_vol_size.first, beampipe_vol_size.second,
                        min_z, max_z, -M_PI, M_PI});
    const auto beampipe_idx = beampipe.index();
    beampipe.set_link(detector_t::sf_finders::id::e_default, 0);

    // This is the beampipe surface
    typename detector_t::surface_type::volume_link_type volume_link{
        beampipe_idx};
    add_cylinder_surface(beampipe_idx, ctx, surfaces, masks, materials,
                         transforms, beampipe_r, min_z, max_z, volume_link,
                         beryllium_tml<scalar>(), 0.8 * unit<scalar>::mm);

    // Get vol sizes in z, including for gap volumes
    std::vector<std::pair<scalar, scalar>> vol_sizes{
        {brl_half_z, edc_lay_sizes[0].first}};
    for (std::size_t i = 0; i < n_edc_layers; ++i) {
        vol_sizes.emplace_back(edc_lay_sizes[i].first, edc_lay_sizes[i].second);
        vol_sizes.emplace_back(edc_lay_sizes[i].second,
                               edc_lay_sizes[i + 1].first);
    }
    vol_sizes.pop_back();

    // negative endcap portals
    dindex link = beampipe_idx;
    for (int i = vol_sizes.size() - 1; i >= 0; --i) {
        volume_link = ++link;
        add_cylinder_surface(beampipe_idx, ctx, surfaces, masks, materials,
                             transforms, edc_inner_r, -vol_sizes[i].second,
                             -vol_sizes[i].first, volume_link, vacuum<scalar>(),
                             0. * unit<scalar>::mm);
    }

    // barrel portals
    volume_link = n_brl_layers <= 0 ? leaving_world : link + 1;
    add_cylinder_surface(beampipe_idx, ctx, surfaces, masks, materials,
                         transforms, edc_inner_r, -brl_half_z, brl_half_z,
                         volume_link, vacuum<scalar>(), 0. * unit<scalar>::mm);

    // positive endcap portals
    link += 7;
    for (std::size_t i = 0; i < vol_sizes.size(); ++i) {
        volume_link = ++link;
        add_cylinder_surface(beampipe_idx, ctx, surfaces, masks, materials,
                             transforms, edc_inner_r, vol_sizes[i].second,
                             vol_sizes[i].first, volume_link, vacuum<scalar>(),
                             0. * unit<scalar>::mm);
    }

    // disc portals
    volume_link = leaving_world;
    add_disc_surface(beampipe_idx, ctx, surfaces, masks, materials, transforms,
                     beampipe_vol_size.first, beampipe_vol_size.second, min_z,
                     volume_link, vacuum<scalar>(), 0. * unit<scalar>::mm);
    add_disc_surface(beampipe_idx, ctx, surfaces, masks, materials, transforms,
                     beampipe_vol_size.first, beampipe_vol_size.second, max_z,
                     volume_link, vacuum<scalar>(), 0. * unit<scalar>::mm);

    det.add_objects_per_volume(ctx, beampipe, surfaces, masks, transforms,
                               materials);
}

/** Helper method for creating a connecting gap volume between endcap and barrel
 *
 * @param det detector the volume should be added to
 * @param resource vecmem memory resource for the temporary containers
 * @param ctx geometry context
 * @param n_brl_layers number of barrel layers that contain modules
 * @param beampipe_idx index of the beampipe volume
 * @param brl_lay_sizes extend of the disc portal surfaces in r
 * @param edc_inner_r inner radius of the gap volume
 * @param edc_outer_r outer radius of the gap volume
 * @param gap_neg_z lower extend of the gap volume in z
 * @param gap_pos_z upper extend of the gap volume in z
 * @param brl_vol_idx index of the first barrel volume (innermost layer)
 * @param edc_vol_idx index of the bordering endcap volume
 */
template <typename detector_t>
inline void add_endcap_barrel_connection(
    detector_t &det, vecmem::memory_resource &resource,
    typename detector_t::geometry_context &ctx, const int side,
    const unsigned int n_brl_layers, const dindex beampipe_idx,
    const std::vector<std::pair<scalar, scalar>> &brl_lay_sizes,
    const scalar edc_inner_r, const scalar edc_outer_r,
    const scalar gap_lower_z, const scalar gap_upper_z, dindex brl_vol_idx,
    const dindex edc_vol_idx) {
    const scalar min_z = std::min(side * gap_lower_z, side * gap_upper_z);
    const scalar max_z = std::max(side * gap_lower_z, side * gap_upper_z);
    scalar edc_disc_z = side < 0 ? min_z : max_z;
    scalar brl_disc_z = side < 0 ? max_z : min_z;

    typename detector_t::surface_container_t surfaces(&resource);
    typename detector_t::mask_container masks(resource);
    typename detector_t::material_container materials(resource);
    typename detector_t::transform_container transforms(resource);

    auto &connector_gap =
        det.new_volume(volume_id::e_cylinder,
                       {edc_inner_r, edc_outer_r, min_z, max_z, -M_PI, M_PI});
    connector_gap.set_link(detector_t::sf_finders::id::e_default, 0);
    dindex connector_gap_idx{det.volumes().back().index()};
    dindex leaving_world = dindex_invalid;

    typename detector_t::surface_type::volume_link_type volume_link = {
        beampipe_idx};
    add_cylinder_surface(connector_gap_idx, ctx, surfaces, masks, materials,
                         transforms, edc_inner_r, min_z, max_z, volume_link,
                         vacuum<scalar>(), 0. * unit<scalar>::mm);
    volume_link = leaving_world;
    add_cylinder_surface(connector_gap_idx, ctx, surfaces, masks, materials,
                         transforms, edc_outer_r, min_z, max_z, volume_link,
                         vacuum<scalar>(), 0. * unit<scalar>::mm);
    volume_link = edc_vol_idx;
    add_disc_surface(connector_gap_idx, ctx, surfaces, masks, materials,
                     transforms, edc_inner_r, edc_outer_r, edc_disc_z,
                     volume_link, vacuum<scalar>(), 0. * unit<scalar>::mm);

    // Get vol sizes in z also for gap volumes
    std::vector<std::pair<scalar, scalar>> vol_sizes;
    for (std::size_t i = 1; i <= n_brl_layers; ++i) {
        vol_sizes.emplace_back(brl_lay_sizes[i].first, brl_lay_sizes[i].second);
        vol_sizes.emplace_back(brl_lay_sizes[i].second,
                               brl_lay_sizes[i + 1].first);
    }

    volume_link = brl_vol_idx;
    for (std::size_t i = 0; i < 2 * n_brl_layers - 1; ++i) {
        volume_link = brl_vol_idx++;
        add_disc_surface(connector_gap_idx, ctx, surfaces, masks, materials,
                         transforms, vol_sizes[i].first, vol_sizes[i].second,
                         brl_disc_z, volume_link, vacuum<scalar>(),
                         0. * unit<scalar>::mm);
    }

    det.add_objects_per_volume(ctx, connector_gap, surfaces, masks, transforms,
                               materials);
}

/** Helper method for creating one of the two endcaps.
 *
 * @param det detector the subdetector should be added to
 * @param resource vecmem memory resource for the temporary containers
 * @param ctx geometry context
 * @param n_layers number of layers that contain modules
 * @param beampipe_idx index of the beampipe outermost volume
 * @param lay_sizes extend of the endcap layers in z direction
 * @param lay_positions position of the endcap layers in z direction
 * @param cfg config struct for module creation
 */
template <typename empty_vol_factory, typename edc_module_factory,
          typename detector_t, typename config_t>
void add_endcap_detector(
    detector_t &det, vecmem::memory_resource &resource,
    typename detector_t::geometry_context &ctx, std::size_t n_layers,
    dindex beampipe_idx,
    const std::vector<std::pair<scalar, scalar>> &lay_sizes,
    const std::vector<scalar> &lay_positions, config_t cfg) {

    // Generate consecutive linking between volumes (all volume_links for every
    // vol.)
    using volume_link_t = typename detector_t::surface_type::volume_link_type;
    std::vector<std::vector<volume_link_t>> volume_links_vec;
    dindex leaving_world = dindex_invalid;
    dindex first_vol_idx = det.volumes().back().index() + 1;
    dindex last_vol_idx = first_vol_idx + 2 * n_layers - 2;
    dindex prev_vol_idx = first_vol_idx - 1;
    dindex next_vol_idx = first_vol_idx + 1;

    for (int i = 0; i < 2 * static_cast<int>(n_layers) - 3; ++i) {
        volume_links_vec.push_back(
            {beampipe_idx, leaving_world, ++prev_vol_idx, ++next_vol_idx});
    }
    // Edge of the world is flipped
    if (cfg.side < 0) {
        volume_links_vec.insert(
            volume_links_vec.begin(),
            {beampipe_idx, leaving_world, leaving_world, first_vol_idx + 1});

        volume_links_vec.push_back(
            {beampipe_idx, leaving_world, last_vol_idx - 1, last_vol_idx + 1});
    } else {
        // For n_layers=1 no extra gap layer is needed
        if (n_layers > 1) {
            volume_links_vec.insert(volume_links_vec.begin(),
                                    {beampipe_idx, leaving_world,
                                     first_vol_idx - 1, first_vol_idx + 1});
        }

        volume_links_vec.push_back(
            {beampipe_idx, leaving_world, last_vol_idx - 1, leaving_world});
    }

    // Get vol sizes in z, including gap volumes
    std::vector<std::pair<scalar, scalar>> vol_sizes{
        {lay_sizes[0].first, lay_sizes[0].second}};
    for (std::size_t i = 1; i < n_layers; ++i) {
        vol_sizes.emplace_back(lay_sizes[i].first, lay_sizes[i - 1].second);
        vol_sizes.emplace_back(lay_sizes[i].first, lay_sizes[i].second);
    }

    edc_module_factory m_factory{cfg};
    empty_vol_factory empty_factory{};

    auto vol_size_itr = vol_sizes.begin();
    auto pos_itr = lay_positions.begin();
    // Reverse iteration for negative endcap
    if (cfg.side < 0) {
        std::advance(vol_size_itr, 2 * n_layers - 2);
        std::advance(pos_itr, n_layers - 1);
    }
    bool is_gap = true;
    for (std::size_t i = 0; i < 2 * n_layers - 1; ++i) {
        // Every second layer is a gap volume
        is_gap = !is_gap;
        if (is_gap) {
            create_cyl_volume(det, resource, ctx, cfg.inner_r, cfg.outer_r,
                              cfg.side * (vol_size_itr + cfg.side * i)->first,
                              cfg.side * (vol_size_itr + cfg.side * i)->second,
                              volume_links_vec[i], empty_factory);

        } else {
            m_factory.cfg.edc_position = *(pos_itr + cfg.side * i / 2);
            create_cyl_volume(det, resource, ctx, cfg.inner_r, cfg.outer_r,
                              cfg.side * (vol_size_itr + cfg.side * i)->first,
                              cfg.side * (vol_size_itr + cfg.side * i)->second,
                              volume_links_vec[i], m_factory);
            add_disc_grid(ctx, det.volumes().back(), det, cfg);
        }
    }
}

/** Helper method for creating the barrel section.
 *
 * @param det detector the subdetector should be added to
 * @param resource vecmem memory resource for the temporary containers
 * @param ctx geometry context
 * @param n_layers number of layers that contain modules
 * @param beampipe_idx index of the beampipe outermost volume
 * @param brl_half_z half length of the barrel section in z direction
 * @param lay_sizes extend of the barrel layers in r direction
 * @param lay_positions position of the barrel layers in r direction
 * @param cfg config struct for module creation
 */
template <typename empty_vol_factory, typename brl_module_factory,
          typename detector_t, typename config_t>
void add_barrel_detector(
    detector_t &det, vecmem::memory_resource &resource,
    typename detector_t::geometry_context &ctx, const unsigned int n_layers,
    dindex beampipe_idx, const scalar brl_half_z,
    const std::vector<std::pair<scalar, scalar>> &lay_sizes,
    const std::vector<scalar> &lay_positions,
    const std::vector<std::pair<int, int>> &m_binning, config_t cfg) {

    // Generate consecutive linking between volumes
    dindex leaving_world = dindex_invalid;
    dindex first_vol_idx = det.volumes().back().index();
    dindex last_vol_idx = first_vol_idx + 2 * n_layers;
    dindex prev_vol_idx = first_vol_idx;
    dindex next_vol_idx = n_layers > 1 ? first_vol_idx + 2 : leaving_world;

    // Leave world, if no endcaps are present
    if (det.volumes().back().index() == 0) {
        first_vol_idx = leaving_world;
        last_vol_idx = leaving_world;
    }

    // First barrel layer is connected to the beampipe
    using volume_link_t = typename detector_t::surface_type::volume_link_type;
    std::vector<std::vector<volume_link_t>> volume_links_vec{
        {beampipe_idx, next_vol_idx, first_vol_idx, last_vol_idx}};

    for (std::size_t i = 1; i < 2 * n_layers - 2; ++i) {
        volume_links_vec.push_back(
            {++prev_vol_idx, ++next_vol_idx, first_vol_idx, last_vol_idx});
    }
    // Last barrel layer leaves detector world
    volume_links_vec.push_back(
        {++prev_vol_idx, leaving_world, first_vol_idx, last_vol_idx});

    // Get vol sizes in z, including gap volumes
    std::vector<std::pair<scalar, scalar>> vol_sizes{
        {lay_sizes[1].first, lay_sizes[1].second}};
    for (std::size_t i = 2; i < n_layers + 1; ++i) {
        vol_sizes.emplace_back(lay_sizes[i].first, lay_sizes[i - 1].second);
        vol_sizes.emplace_back(lay_sizes[i].first, lay_sizes[i].second);
    }

    brl_module_factory m_factory{cfg};
    empty_vol_factory empty_factory{};
    bool is_gap = true;
    for (unsigned int i = 0; i < 2 * n_layers - 1; ++i) {
        unsigned int j = (i + 2) / 2;
        // Every second layer is a gap volume
        is_gap = !is_gap;
        if (is_gap) {
            create_cyl_volume(det, resource, ctx, vol_sizes[i].first,
                              vol_sizes[i].second, -brl_half_z, brl_half_z,
                              volume_links_vec[i], empty_factory);
        } else {
            m_factory.cfg.m_binning = m_binning[j];
            m_factory.cfg.layer_r = lay_positions[j];
            create_cyl_volume(det, resource, ctx, vol_sizes[i].first,
                              vol_sizes[i].second, -brl_half_z, brl_half_z,
                              volume_links_vec[i], m_factory);
            // add_cylinder_grid(ctx, det.volumes().back(), det, cfg);
        }
    }
}

}  // namespace

/** Builds a detray geometry that contains the innermost tml layers. The number
 *  of barrel and endcap layers can be chosen, but all barrel layers should be
 *  present when an endcap detector is built to have the barrel region radius
 *  match the endcap diameter.
 *
 * @param n_brl_layers number of pixel barrel layer to build (max 4)
 * @param n_edc_layers number of pixel endcap discs to build (max 7)
 *
 * @returns a complete detector object
 */
template <typename container_t = host_container_types>
auto create_toy_geometry(
    vecmem::memory_resource &resource,
    covfie::field<detector_registry::toy_detector::bfield_backend_t> &&bfield,
    std::size_t n_brl_layers = 4, std::size_t n_edc_layers = 3) {

    // detector type
    using detector_t =
        detector<detector_registry::toy_detector, covfie::field, container_t>;

    /// Leaving world
    constexpr dindex leaving_world{dindex_invalid};

    //
    // barrel
    //
    constexpr scalar brl_half_z{500.};
    const std::vector<scalar> brl_positions = {19., 32., 72., 116., 172.};
    const std::vector<std::pair<scalar, scalar>> brl_lay_sizes = {
        {0., 27.}, {27., 38.}, {64., 80.}, {108., 124.}, {164., 180.}};
    const std::vector<std::pair<int, int>> brl_binning = {
        {0., 0.}, {16, 14}, {32, 14}, {52, 14}, {78, 14}};
    // module parameters
    struct brl_m_config {
        scalar m_half_x{8.4};
        scalar m_half_y{36.};
        scalar m_tilt_phi{0.14};  // 0.145;
        scalar layer_r{32.};
        scalar m_radial_stagger{0.5};  // 2.;
        scalar m_long_overlap{2.};     // 5.;
        std::pair<std::size_t, std::size_t> m_binning = {16, 14};
        material<scalar> mat = silicon_tml<scalar>();
        scalar thickness = 0.15 * unit<scalar>::mm;
    };

    //
    // endcaps
    //
    const std::vector<scalar> edc_positions = {600.,  700.,  820., 960.,
                                               1100., 1300., 1500.};
    const std::vector<std::pair<scalar, scalar>> edc_lay_sizes = {
        {595., 605.},   {695., 705.},   {815., 825.},  {955., 965.},
        {1095., 1105.}, {1295., 1305.}, {1495., 1505.}};
    // module params
    struct edc_m_config {
        int side{1};
        scalar inner_r{27.};
        scalar outer_r{180.};
        scalar edc_position{600.};
        scalar ring_stagger{1.0};
        // Parameters for both rings of modules
        std::vector<scalar> m_phi_stagger = {4.0, 4.0};
        std::vector<scalar> m_phi_sub_stagger = {0.5, 0.5};
        std::vector<std::size_t> disc_binning = {40, 68};
        std::vector<scalar> m_half_y = {36., 36.};
        // std::vector<scalar> m_half_x_min_y = {8.4, 8.4};
        // std::vector<scalar> m_half_x_max_y = {10.1, 10.1};
        std::vector<scalar> m_half_x_min_y = {8.4, 8.4};
        std::vector<scalar> m_half_x_max_y = {12.4, 12.4};
        std::vector<scalar> m_tilt = {0., 0.};
        material<scalar> mat = silicon_tml<scalar>();
        scalar thickness = 0.15 * unit<scalar>::mm;
    };

    // Don't create modules in gap volume
    struct empty_vol_factory {
        void operator()(
            typename detector_t::geometry_context & /*ctx*/,
            typename detector_t::volume_type & /*volume*/,
            typename detector_t::surface_container_t & /*surfaces*/,
            typename detector_t::mask_container & /*masks*/,
            typename detector_t::material_container & /*materials*/,
            typename detector_t::transform_container & /*transforms*/) {}
    };

    // Fills volume with barrel layer
    struct brl_module_factory {
        brl_m_config cfg;

        void operator()(typename detector_t::geometry_context &ctx,
                        typename detector_t::volume_type &volume,
                        typename detector_t::surface_container_t &surfaces,
                        typename detector_t::mask_container &masks,
                        typename detector_t::material_container &materials,
                        typename detector_t::transform_container &transforms) {
            create_barrel_modules(ctx, volume, surfaces, masks, materials,
                                  transforms, cfg);
        }
    };

    // Fills volume with endcap rings
    struct edc_module_factory {
        edc_m_config cfg;

        void operator()(typename detector_t::geometry_context &ctx,
                        typename detector_t::volume_type &volume,
                        typename detector_t::surface_container_t &surfaces,
                        typename detector_t::mask_container &masks,
                        typename detector_t::material_container &materials,
                        typename detector_t::transform_container &transforms) {
            create_endcap_modules(ctx, volume, surfaces, masks, materials,
                                  transforms, cfg);
        }
    };

    // create empty detector
    detector_t det(resource, std::move(bfield));

    // geometry context object
    typename detector_t::geometry_context ctx0{};

    brl_m_config brl_config{};
    edc_m_config edc_config{};

    if (n_edc_layers > edc_positions.size()) {
        throw std::invalid_argument(
            "ERROR: Too many endcap layers requested (max " +
            std::to_string(edc_positions.size()) + ")!");
    }
    if (n_brl_layers > brl_positions.size() - 1) {
        throw std::invalid_argument(
            "ERROR: Too many barrel layers requested (max " +
            std::to_string(brl_positions.size() - 1) + ")!");
    }
    // the radius of the endcaps and  the barrel section need to match
    if (n_edc_layers > 0 and
        std::fabs(brl_lay_sizes[n_brl_layers].second - edc_config.outer_r) >
            std::numeric_limits<scalar>::epsilon()) {
        throw std::invalid_argument(
            "ERROR: Barrel and endcap radii do not match!");
    }

    // beampipe
    dindex beampipe_idx = 0;
    add_beampipe(det, resource, ctx0, n_edc_layers, n_brl_layers, edc_lay_sizes,
                 brl_lay_sizes[0], brl_positions[0], brl_half_z,
                 edc_config.inner_r);

    if (n_edc_layers > 0) {
        edc_config.side = -1;
        // negative endcap layers
        add_endcap_detector<empty_vol_factory, edc_module_factory>(
            det, resource, ctx0, n_edc_layers, beampipe_idx, edc_lay_sizes,
            edc_positions, edc_config);

        // gap volume that connects barrel and neg. endcap
        dindex prev_vol_idx = det.volumes().back().index();
        prev_vol_idx = prev_vol_idx == 0 ? leaving_world : prev_vol_idx;
        dindex next_vol_idx = n_brl_layers == 0
                                  ? leaving_world
                                  : det.volumes().back().index() + 2;

        add_endcap_barrel_connection(
            det, resource, ctx0, edc_config.side, n_brl_layers, beampipe_idx,
            brl_lay_sizes, edc_config.inner_r, edc_config.outer_r,
            edc_lay_sizes[0].first, brl_half_z, next_vol_idx, prev_vol_idx);
    }
    if (n_brl_layers > 0) {
        // barrel
        add_barrel_detector<empty_vol_factory, brl_module_factory>(
            det, resource, ctx0, n_brl_layers, beampipe_idx, brl_half_z,
            brl_lay_sizes, brl_positions, brl_binning, brl_config);
    }
    if (n_edc_layers > 0) {
        // gap layer that connects barrel and pos. endcap
        edc_config.side = 1.;
        // innermost barrel layer volume id
        dindex prev_vol_idx =
            n_brl_layers == 0 ? leaving_world : 2 * n_edc_layers + 1;
        dindex next_vol_idx = prev_vol_idx == 1
                                  ? leaving_world
                                  : det.volumes().back().index() + 2;

        add_endcap_barrel_connection(
            det, resource, ctx0, edc_config.side, n_brl_layers, beampipe_idx,
            brl_lay_sizes, edc_config.inner_r, edc_config.outer_r, brl_half_z,
            edc_lay_sizes[0].first, prev_vol_idx, next_vol_idx);

        // positive endcap layers
        add_endcap_detector<empty_vol_factory, edc_module_factory>(
            det, resource, ctx0, n_edc_layers, beampipe_idx, edc_lay_sizes,
            edc_positions, edc_config);
    }

    return det;
}

/** Wrapper for create_toy_geometry with constant zero bfield.
 */
template <typename container_t = host_container_types>
auto create_toy_geometry(vecmem::memory_resource &resource,
                         std::size_t n_brl_layers = 4,
                         std::size_t n_edc_layers = 3) {
    return create_toy_geometry<container_t>(
        resource,
        covfie::field<detector_registry::toy_detector::bfield_backend_t>{
            detector_registry::toy_detector::bfield_backend_t::configuration_t{
                0.f, 0.f, 0.f}},
        n_brl_layers, n_edc_layers);
}

}  // namespace detray