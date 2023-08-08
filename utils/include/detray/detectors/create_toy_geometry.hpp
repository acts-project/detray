/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/detector_helper.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/tools/grid_builder.hpp"
#include "detray/tools/volume_builder.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/field.hpp>

// System include(s)
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace detray {

namespace {

using transform3_t = __plugin::transform3<detray::scalar>;
using point3 = __plugin::point3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;
using point2 = __plugin::point2<detray::scalar>;

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
    using nav_link_t = typename surface_type::navigation_link;
    using mask_id = typename surface_type::mask_id;
    using mask_link_type = typename surface_type::mask_link;
    using material_id = typename surface_type::material_id;
    using material_link_type = typename surface_type::material_link;

    constexpr auto rectangle_id = mask_id::e_rectangle2;
    constexpr auto slab_id = material_id::e_slab;

    auto volume_idx = vol.index();
    auto mask_volume_link{static_cast<nav_link_t>(volume_idx)};

    // Create the module centers

    // surface grid bins
    unsigned int n_phi_bins = cfg.m_binning.first;
    unsigned int n_z_bins = cfg.m_binning.second;
    // module positions
    detray::dvector<point3> m_centers;
    m_centers.reserve(n_phi_bins * n_z_bins);

    // prep work
    scalar phi_step{2.f * constant<scalar>::pi /
                    static_cast<scalar>(n_phi_bins)};
    scalar min_phi{-constant<scalar>::pi + 0.5f * phi_step};

    scalar z_start{-0.5f * static_cast<scalar>(n_z_bins - 1u) *
                   (2.f * cfg.m_half_y - cfg.m_long_overlap)};
    scalar z_step{(std::abs(z_start) - z_start) /
                  static_cast<scalar>(n_z_bins - 1)};

    // loop over the z bins
    for (unsigned int z_bin = 0u; z_bin < n_z_bins; ++z_bin) {
        // prepare z and r
        scalar m_z{z_start + static_cast<scalar>(z_bin) * z_step};
        scalar m_r{(z_bin % 2u) != 0u
                       ? cfg.layer_r - 0.5f * cfg.m_radial_stagger
                       : cfg.layer_r + 0.5f * cfg.m_radial_stagger};
        for (unsigned int phiBin = 0u; phiBin < n_phi_bins; ++phiBin) {
            // calculate the current phi value
            scalar m_phi{min_phi + static_cast<scalar>(phiBin) * phi_step};
            m_centers.push_back(
                point3{m_r * math_ns::cos(m_phi), m_r * std::sin(m_phi), m_z});
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
        vector3 m_local_z{math_ns::cos(m_phi + cfg.m_tilt_phi),
                          std::sin(m_phi + cfg.m_tilt_phi), 0.f};
        // Local x axis the normal to local y,z
        vector3 m_local_x{-std::sin(m_phi + cfg.m_tilt_phi),
                          math_ns::cos(m_phi + cfg.m_tilt_phi), 0.f};

        // Create the module transform
        transforms.emplace_back(ctx, m_center, m_local_z, m_local_x);
    }
}

/// Helper function that creates a surface grid of rectangular barrel modules.
///
/// @param surfaces_grid the grid to be created with proper axes
/// @param resource vecmem memory resource
/// @param cfg config struct for module creation
template <
    typename detector_t, typename config_t, typename factory_t,
    std::enable_if_t<
        std::is_invocable_v<factory_t, typename detector_t::geometry_context &,
                            typename detector_t::volume_type &,
                            typename detector_t::surface_container_t &,
                            typename detector_t::mask_container &,
                            typename detector_t::material_container &,
                            typename detector_t::transform_container &>,
        bool> = true>
inline void add_cylinder_grid(const typename detector_t::geometry_context &ctx,
                              vecmem::memory_resource &resource,
                              typename detector_t::volume_type &vol,
                              detector_t &det, const config_t &cfg,
                              factory_t &module_factory) {
    // Get relevant ids
    using geo_obj_ids = typename detector_t::geo_obj_ids;

    constexpr auto cyl_id = detector_t::masks::id::e_portal_cylinder2;
    constexpr auto grid_id = detector_t::sf_finders::id::e_cylinder2_grid;

    using cyl_grid_t =
        typename detector_t::surface_container::template get_type<grid_id>;
    auto gbuilder =
        grid_builder<detector_t, cyl_grid_t, detray::detail::fill_by_pos>{};

    // The disc portals are at the end of the portal range by construction
    auto portal_mask_idx = (det.portals(vol).end() - 3u)->mask().index();
    const auto &cyl_mask =
        det.mask_store().template get<cyl_id>().at(portal_mask_idx);

    // Create the sensitive surfaces
    typename detector_t::surface_container_t surfaces(&resource);
    typename detector_t::mask_container masks(resource);
    typename detector_t::material_container materials(resource);
    typename detector_t::transform_container transforms(resource);
    module_factory(ctx, vol, surfaces, masks, materials, transforms);

    // Iterate the surfaces and update their links
    const auto trf_offset{det.transform_store().size(ctx)};
    auto sf_offset{det.n_surfaces()};
    for (auto &sf_desc : surfaces) {
        // Make sure the volume was constructed correctly
        assert(sf_desc.volume() < det.volumes().size());

        // Update the surface links accroding to number of data in detector
        const auto sf = surface{det, sf_desc};
        sf.template visit_mask<detail::mask_index_update>(sf_desc);
        sf.template visit_material<detail::material_index_update>(sf_desc);
        sf_desc.update_transform(trf_offset);
        sf_desc.set_index(sf_offset++);

        // Copy surface descriptor into global lookup
        det.add_surface_to_lookup(sf_desc);
    }

    // Add new grid to the detector
    gbuilder.init_grid(cyl_mask, {cfg.m_binning.first, cfg.m_binning.second});
    gbuilder.fill_grid(detector_volume{det, vol}, surfaces, transforms, masks,
                       ctx);
    assert(gbuilder.get().all().size() == surfaces.size());

    // Add transforms, masks and material to detector
    det.append_masks(std::move(masks));
    det.append_transforms(std::move(transforms));
    det.append_materials(std::move(materials));
    det.surface_store().template push_back<grid_id>(gbuilder.get());
    vol.template set_link<geo_obj_ids::e_sensitive>(
        grid_id, det.surface_store().template size<grid_id>() - 1u);
}

/// Helper function that creates a surface grid of trapezoidal endcap modules.
///
/// @param vol the detector volume that should be equipped with a grid
/// @param surfaces_grid the grid to be created with proper axes
/// @param resource vecmem memory resource
/// @param cfg config struct for module creation
template <
    typename detector_t, typename config_t, typename factory_t,
    std::enable_if_t<
        std::is_invocable_v<factory_t, typename detector_t::geometry_context &,
                            typename detector_t::volume_type &,
                            typename detector_t::surface_container_t &,
                            typename detector_t::mask_container &,
                            typename detector_t::material_container &,
                            typename detector_t::transform_container &>,
        bool> = true>
inline void add_disc_grid(const typename detector_t::geometry_context &ctx,
                          vecmem::memory_resource &resource,
                          typename detector_t::volume_type &vol,
                          detector_t &det, const config_t &cfg,
                          factory_t &module_factory) {
    // Get relevant ids
    using geo_obj_ids = typename detector_t::geo_obj_ids;

    constexpr auto disc_id = detector_t::masks::id::e_portal_ring2;
    constexpr auto grid_id = detector_t::sf_finders::id::e_disc_grid;

    using disc_grid_t =
        typename detector_t::surface_container::template get_type<grid_id>;
    auto gbuilder =
        grid_builder<detector_t, disc_grid_t, detray::detail::fill_by_pos>{};

    // The disc portals are at the end of the portal range by construction
    auto portal_mask_idx = det.portals(vol).back().mask().index();
    const auto &disc_mask =
        det.mask_store().template get<disc_id>().at(portal_mask_idx);

    // Create the sensitive surfaces
    typename detector_t::surface_container_t surfaces(&resource);
    typename detector_t::mask_container masks(resource);
    typename detector_t::material_container materials(resource);
    typename detector_t::transform_container transforms(resource);
    module_factory(ctx, vol, surfaces, masks, materials, transforms);

    // Iterate the surfaces and update their links
    const auto trf_offset{det.transform_store().size(ctx)};
    auto sf_offset{det.n_surfaces()};
    for (auto &sf_desc : surfaces) {
        // Make sure the volume was constructed correctly
        assert(sf_desc.volume() < det.volumes().size());

        // Update the surface links accroding to number of data in detector
        const auto sf = surface{det, sf_desc};
        sf.template visit_mask<detail::mask_index_update>(sf_desc);
        sf.template visit_material<detail::material_index_update>(sf_desc);
        sf_desc.update_transform(trf_offset);
        sf_desc.set_index(sf_offset++);

        // Copy surface descriptor into global lookup
        det.add_surface_to_lookup(sf_desc);
    }

    // Add new grid to the detector
    gbuilder.init_grid(disc_mask,
                       {cfg.disc_binning.size(), cfg.disc_binning.front()});
    gbuilder.fill_grid(detector_volume{det, vol}, surfaces, transforms, masks,
                       ctx);
    assert(gbuilder.get().all().size() == surfaces.size());

    // Add transforms, masks and material to detector
    det.append_masks(std::move(masks));
    det.append_transforms(std::move(transforms));
    det.append_materials(std::move(materials));
    det.surface_store().template push_back<grid_id>(gbuilder.get());
    vol.template set_link<geo_obj_ids::e_sensitive>(
        grid_id, det.surface_store().template size<grid_id>() - 1u);
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
                                  scalar phi_sub_stagger,
                                  unsigned int n_phi_bins) {
    // create and fill the positions
    std::vector<vector3> r_positions;
    r_positions.reserve(n_phi_bins);

    // prep work
    scalar phi_step{2.f * constant<scalar>::pi /
                    static_cast<scalar>(n_phi_bins)};
    scalar min_phi{-constant<scalar>::pi + 0.5f * phi_step};

    for (unsigned int iphi = 0u; iphi < n_phi_bins; ++iphi) {
        // if we have a phi sub stagger presents
        scalar rzs{0.f};
        // phi stagger affects 0 vs 1, 2 vs 3 ... etc
        // -> only works if it is a %4
        // phi sub stagger affects 2 vs 4, 1 vs 3 etc.
        if (phi_sub_stagger != 0.f && !(n_phi_bins % 4u)) {
            // switch sides
            if (!(iphi % 4u)) {
                rzs = phi_sub_stagger;
            } else if (!((iphi + 1u) % 4u)) {
                rzs = -phi_sub_stagger;
            }
        }
        // the module phi
        scalar phi{min_phi + static_cast<scalar>(iphi) * phi_step};
        // main z position depending on phi bin
        scalar rz{iphi % 2u ? z - 0.5f * phi_stagger : z + 0.5f * phi_stagger};
        r_positions.push_back(vector3{radius * math_ns::cos(phi),
                                      radius * std::sin(phi), rz + rzs});
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
    using nav_link_t = typename surface_type::navigation_link;
    using mask_id = typename surface_type::mask_id;
    using mask_link_type = typename surface_type::mask_link;
    using material_id = typename surface_type::material_id;
    using material_link_type = typename surface_type::material_link;

    constexpr auto trapezoid_id = mask_id::e_trapezoid2;
    constexpr auto slab_id = material_id::e_slab;

    auto volume_idx = vol.index();
    auto mask_volume_link{static_cast<nav_link_t>(volume_idx)};

    // calculate the radii of the rings
    std::vector<scalar> radii;
    // calculate the radial borders
    // std::vector<scalar> radial_boarders;
    // the radial span of the disc
    scalar delta_r{cfg.outer_r - cfg.inner_r};

    // Only one ring
    if (cfg.disc_binning.size() == 1u) {
        radii.push_back(0.5f * (cfg.inner_r + cfg.outer_r));
        // radial_boarders = {inner_r, outer_r};
    } else {
        // sum up the total length of the modules along r
        scalar tot_length{0.f};
        for (auto &m_hlength : cfg.m_half_y) {
            tot_length += 2.f * m_hlength + 0.5f;
        }
        // now calculate the overlap (equal pay)
        scalar r_overlap{(tot_length - delta_r) /
                         static_cast<scalar>(cfg.m_half_y.size() - 1u)};
        // and now fill the radii and gaps
        scalar prev_r{cfg.inner_r};
        scalar prev_hl{0.f};
        scalar prev_ol{0.f};
        // remember the radial boarders
        // radial_boarders.push_back(inner_r);
        for (auto &m_hlength : cfg.m_half_y) {
            // calculate the radius
            radii.push_back(prev_r + prev_hl - prev_ol + m_hlength);
            prev_r = radii.back();
            prev_ol = r_overlap;
            prev_hl = m_hlength;
            // and register the radial boarder
            // radial_boarders.push_back(prev_r + 2.f * prev_hl -
            // 0.5f * prev_ol);
        }
    }

    // now build the modules in every ring
    for (unsigned int ir = 0u; ir < radii.size(); ++ir) {
        // generate the z value
        // convention inner ring is closer to origin : makes sense
        scalar rz{radii.size() == 1u
                      ? cfg.edc_position
                      : (ir % 2u ? cfg.edc_position + 0.5f * cfg.ring_stagger
                                 : cfg.edc_position - 0.5f * cfg.ring_stagger)};
        // fill the ring module positions
        scalar ps_stagger{
            cfg.m_phi_sub_stagger.size() ? cfg.m_phi_sub_stagger[ir] : 0.f};

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
                1.f / (2.f * cfg.m_half_y[ir]));

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
            vector3 m_local_y{math_ns::cos(m_phi), std::sin(m_phi), 0.f};
            // take different axis to have the same readout direction
            vector3 m_local_z{0.f, 0.f, static_cast<scalar>(cfg.side)};
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
    typename detector_t::geometry_context &ctx,
    typename detector_t::name_map &names, const unsigned int n_edc_layers,
    const unsigned int n_brl_layers,
    const std::vector<std::pair<scalar, scalar>> &edc_lay_sizes,
    const std::pair<scalar, scalar> &beampipe_vol_size, const scalar beampipe_r,
    const scalar brl_half_z, const scalar edc_inner_r) {

    // Get the object types handled by the volume
    using object_id = typename detector_t::volume_type::object_id;
    using nav_link_t = typename detector_t::surface_type::navigation_link;
    constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};
    constexpr auto cyl_id = detector_t::masks::id::e_portal_cylinder2;

    scalar max_z{n_edc_layers == 0u ? brl_half_z
                                    : edc_lay_sizes[n_edc_layers - 1u].second};
    scalar min_z{-max_z};

    typename detector_t::surface_container_t surfaces(&resource);
    typename detector_t::mask_container masks(resource);
    typename detector_t::material_container materials(resource);
    typename detector_t::transform_container transforms(resource);

    auto &beampipe = det.new_volume(volume_id::e_cylinder);
    const auto beampipe_idx = beampipe.index();
    names[beampipe_idx + 1u] = "beampipe_" + std::to_string(beampipe_idx);
    beampipe.set_transform(det.transform_store().size());
    det.transform_store().emplace_back(ctx);

    beampipe.template set_link<object_id::e_portal>(
        detector_t::sf_finders::id::e_default, 0u);
    beampipe.template set_link<object_id::e_passive>(
        detector_t::sf_finders::id::e_default, 0u);

    // This is the beampipe surface
    dindex volume_link{beampipe_idx};

    // Get vol sizes in z, including for gap volumes
    std::vector<std::pair<scalar, scalar>> vol_sizes{
        {brl_half_z, edc_lay_sizes[0].first}};
    for (unsigned int i = 0u; i < n_edc_layers; ++i) {
        vol_sizes.emplace_back(edc_lay_sizes[i].first, edc_lay_sizes[i].second);
        vol_sizes.emplace_back(edc_lay_sizes[i].second,
                               edc_lay_sizes[i + 1].first);
    }
    vol_sizes.pop_back();

    // negative endcap portals
    dindex link = beampipe_idx;
    for (int i = static_cast<int>(vol_sizes.size()) - 1; i >= 0; --i) {
        volume_link = ++link;
        detail::detector_helper<transform3_t>().add_cylinder_surface<cyl_id>(
            beampipe_idx, ctx, surfaces, masks, materials, transforms,
            edc_inner_r, -vol_sizes[static_cast<unsigned int>(i)].second,
            -vol_sizes[static_cast<unsigned int>(i)].first, volume_link,
            vacuum<scalar>(), 0.f * unit<scalar>::mm);
    }

    // barrel portals
    volume_link = n_brl_layers <= 0u ? leaving_world : link + 1u;
    detail::detector_helper<transform3_t>().add_cylinder_surface<cyl_id>(
        beampipe_idx, ctx, surfaces, masks, materials, transforms, edc_inner_r,
        -brl_half_z, brl_half_z, volume_link, vacuum<scalar>(),
        0.f * unit<scalar>::mm);

    // positive endcap portals
    link += 7u;
    for (unsigned int i = 0u; i < vol_sizes.size(); ++i) {
        volume_link = ++link;
        detail::detector_helper<transform3_t>().add_cylinder_surface<cyl_id>(
            beampipe_idx, ctx, surfaces, masks, materials, transforms,
            edc_inner_r, vol_sizes[i].second, vol_sizes[i].first, volume_link,
            vacuum<scalar>(), 0.f * unit<scalar>::mm);
    }

    // disc portals
    volume_link = leaving_world;
    detail::detector_helper<transform3_t>().add_disc_surface(
        beampipe_idx, ctx, surfaces, masks, materials, transforms,
        beampipe_vol_size.first, beampipe_vol_size.second, min_z, volume_link,
        vacuum<scalar>(), 0.f * unit<scalar>::mm);
    detail::detector_helper<transform3_t>().add_disc_surface(
        beampipe_idx, ctx, surfaces, masks, materials, transforms,
        beampipe_vol_size.first, beampipe_vol_size.second, max_z, volume_link,
        vacuum<scalar>(), 0.f * unit<scalar>::mm);

    // This is the beampipe surface
    detail::detector_helper<transform3_t>()
        .add_cylinder_surface<detector_t::masks::id::e_cylinder2>(
            beampipe_idx, ctx, surfaces, masks, materials, transforms,
            beampipe_r, min_z, max_z, beampipe_idx, beryllium_tml<scalar>(),
            0.8f * unit<scalar>::mm);

    det.add_objects_per_volume(ctx, beampipe, surfaces, masks, transforms,
                               materials);
}

/** Helper method for creating a connecting gap volume between endcap and barrel
 *
 * @param det detector the volume should be added to
 * @param resource vecmem memory resource for the temporary containers
 * @param ctx geometry context
 * @param side positive vs negative endcap
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
    typename detector_t::geometry_context &ctx,
    typename detector_t::name_map &names, const int side,
    const unsigned int n_brl_layers, const dindex beampipe_idx,
    const std::vector<std::pair<scalar, scalar>> &brl_lay_sizes,
    const scalar edc_inner_r, const scalar edc_outer_r,
    const scalar gap_lower_z, const scalar gap_upper_z, dindex brl_vol_idx,
    const dindex edc_vol_idx) {

    using nav_link_t = typename detector_t::surface_type::navigation_link;
    constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};
    constexpr auto cyl_id = detector_t::masks::id::e_portal_cylinder2;

    const scalar sign{static_cast<scalar>(side)};
    const scalar min_z{std::min(sign * gap_lower_z, sign * gap_upper_z)};
    const scalar max_z{std::max(sign * gap_lower_z, sign * gap_upper_z)};
    const scalar edc_disc_z{side < 0 ? min_z : max_z};
    const scalar brl_disc_z{side < 0 ? max_z : min_z};

    typename detector_t::surface_container_t surfaces(&resource);
    typename detector_t::mask_container masks(resource);
    typename detector_t::material_container materials(resource);
    typename detector_t::transform_container transforms(resource);

    auto &connector_gap = det.new_volume(volume_id::e_cylinder,
                                         {detector_t::sf_finders::id::e_default,
                                          detail::invalid_value<dindex>()});
    dindex connector_gap_idx{det.volumes().back().index()};
    names[connector_gap_idx + 1u] =
        "connector_gap_" + std::to_string(connector_gap_idx);
    connector_gap.set_transform(det.transform_store().size());
    // translation of the connector gap
    point3 t{0.f, 0.f, 0.5f * (max_z + min_z)};
    det.transform_store().emplace_back(ctx, t);

    dindex volume_link{beampipe_idx};
    detail::detector_helper<transform3_t>().add_cylinder_surface<cyl_id>(
        connector_gap_idx, ctx, surfaces, masks, materials, transforms,
        edc_inner_r, min_z, max_z, volume_link, vacuum<scalar>(),
        0.f * unit<scalar>::mm);
    volume_link = leaving_world;
    detail::detector_helper<transform3_t>().add_cylinder_surface<cyl_id>(
        connector_gap_idx, ctx, surfaces, masks, materials, transforms,
        edc_outer_r, min_z, max_z, volume_link, vacuum<scalar>(),
        0.f * unit<scalar>::mm);
    volume_link = edc_vol_idx;
    detail::detector_helper<transform3_t>().add_disc_surface(
        connector_gap_idx, ctx, surfaces, masks, materials, transforms,
        edc_inner_r, edc_outer_r, edc_disc_z, volume_link, vacuum<scalar>(),
        0.f * unit<scalar>::mm);

    // Get vol sizes in z also for gap volumes
    std::vector<std::pair<scalar, scalar>> vol_sizes;
    for (unsigned int i = 1u; i <= n_brl_layers; ++i) {
        vol_sizes.emplace_back(brl_lay_sizes[i].first, brl_lay_sizes[i].second);
        vol_sizes.emplace_back(brl_lay_sizes[i].second,
                               brl_lay_sizes[i + 1u].first);
    }

    volume_link = brl_vol_idx;
    for (unsigned int i = 0u; i < 2u * n_brl_layers - 1u; ++i) {
        volume_link = brl_vol_idx++;
        detail::detector_helper<transform3_t>().add_disc_surface(
            connector_gap_idx, ctx, surfaces, masks, materials, transforms,
            vol_sizes[i].first, vol_sizes[i].second, brl_disc_z, volume_link,
            vacuum<scalar>(), 0.f * unit<scalar>::mm);
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
template <typename edc_module_factory, typename detector_t, typename config_t>
void add_endcap_detector(
    detector_t &det, vecmem::memory_resource &resource,
    typename detector_t::geometry_context &ctx,
    typename detector_t::name_map &names, dindex n_layers, dindex beampipe_idx,
    const std::vector<std::pair<scalar, scalar>> &lay_sizes,
    const std::vector<scalar> &lay_positions, config_t cfg) {
    using nav_link_t = typename detector_t::surface_type::navigation_link;
    constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};

    // Generate consecutive linking between volumes (all volume_links for every
    // vol.)
    std::vector<std::vector<dindex>> volume_links_vec;
    dindex first_vol_idx = det.volumes().back().index() + 1u;
    dindex last_vol_idx = first_vol_idx + 2u * n_layers - 2u;
    dindex prev_vol_idx = first_vol_idx - 1u;
    dindex next_vol_idx = first_vol_idx + 1u;

    for (int i = 0; i < 2 * static_cast<int>(n_layers) - 3; ++i) {
        volume_links_vec.push_back(
            {beampipe_idx, leaving_world, ++prev_vol_idx, ++next_vol_idx});
    }
    // Edge of the world is flipped
    if (cfg.side < 0) {
        volume_links_vec.insert(
            volume_links_vec.begin(),
            {beampipe_idx, leaving_world, leaving_world, first_vol_idx + 1u});

        volume_links_vec.push_back({beampipe_idx, leaving_world,
                                    last_vol_idx - 1u, last_vol_idx + 1u});
    } else {
        // For n_layers=1 no extra gap layer is needed
        if (n_layers > 1u) {
            volume_links_vec.insert(volume_links_vec.begin(),
                                    {beampipe_idx, leaving_world,
                                     first_vol_idx - 1u, first_vol_idx + 1u});
        }

        volume_links_vec.push_back(
            {beampipe_idx, leaving_world, last_vol_idx - 1u, leaving_world});
    }

    // Get vol sizes in z, including gap volumes
    std::vector<std::pair<scalar, scalar>> vol_sizes{
        {lay_sizes[0].first, lay_sizes[0].second}};
    for (unsigned int i = 1u; i < n_layers; ++i) {
        vol_sizes.emplace_back(lay_sizes[i].first, lay_sizes[i - 1u].second);
        vol_sizes.emplace_back(lay_sizes[i].first, lay_sizes[i].second);
    }

    edc_module_factory m_factory{cfg};

    auto vol_size_itr = vol_sizes.begin();
    auto pos_itr = lay_positions.begin();
    // Reverse iteration for negative endcap
    if (cfg.side < 0) {
        std::advance(vol_size_itr, 2u * n_layers - 2u);
        std::advance(pos_itr, n_layers - 1u);
    }
    bool is_gap = true;
    for (int i = 0; i < 2 * static_cast<int>(n_layers) - 1; ++i) {
        // Every second layer is a gap volume
        is_gap = !is_gap;
        const scalar sign{static_cast<scalar>(cfg.side)};
        if (is_gap) {
            detail::detector_helper<transform3_t>().create_cyl_volume(
                det, resource, ctx, cfg.inner_r, cfg.outer_r,
                sign * (vol_size_itr + cfg.side * i)->first,
                sign * (vol_size_itr + cfg.side * i)->second,
                volume_links_vec[static_cast<unsigned int>(i)]);

            dindex vol_idx = det.volumes().back().index();
            names[vol_idx + 1u] = "gap_" + std::to_string(vol_idx);

        } else {
            m_factory.cfg.edc_position = *(pos_itr + cfg.side * i / 2);
            detail::detector_helper<transform3_t>().create_cyl_volume(
                det, resource, ctx, cfg.inner_r, cfg.outer_r,
                sign * (vol_size_itr + cfg.side * i)->first,
                sign * (vol_size_itr + cfg.side * i)->second,
                volume_links_vec[static_cast<unsigned int>(i)]);

            dindex vol_idx = det.volumes().back().index();
            names[vol_idx + 1u] = "endcap_" + std::to_string(vol_idx);

            add_disc_grid(ctx, resource, det.volumes().back(), det, cfg,
                          m_factory);
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
template <typename brl_module_factory, typename detector_t, typename config_t>
void add_barrel_detector(
    detector_t &det, vecmem::memory_resource &resource,
    typename detector_t::geometry_context &ctx,
    typename detector_t::name_map &names, const unsigned int n_layers,
    dindex beampipe_idx, const scalar brl_half_z,
    const std::vector<std::pair<scalar, scalar>> &lay_sizes,
    const std::vector<scalar> &lay_positions,
    const std::vector<std::pair<int, int>> &m_binning, config_t cfg) {

    using nav_link_t = typename detector_t::surface_type::navigation_link;
    constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};

    // Generate consecutive linking between volumes
    dindex first_vol_idx = det.volumes().back().index();
    dindex last_vol_idx = first_vol_idx + 2u * n_layers;
    dindex prev_vol_idx = first_vol_idx;
    dindex next_vol_idx = n_layers > 1u ? first_vol_idx + 2u : leaving_world;

    // Leave world, if no endcaps are present
    if (det.volumes().back().index() == 0u) {
        first_vol_idx = leaving_world;
        last_vol_idx = leaving_world;
    }

    // First barrel layer is connected to the beampipe
    std::vector<std::vector<dindex>> volume_links_vec{
        {beampipe_idx, next_vol_idx, first_vol_idx, last_vol_idx}};

    for (unsigned int i = 1u; i < 2u * n_layers - 2u; ++i) {
        volume_links_vec.push_back(
            {++prev_vol_idx, ++next_vol_idx, first_vol_idx, last_vol_idx});
    }
    // Last barrel layer leaves detector world
    volume_links_vec.push_back(
        {++prev_vol_idx, leaving_world, first_vol_idx, last_vol_idx});

    // Get vol sizes in z, including gap volumes
    std::vector<std::pair<scalar, scalar>> vol_sizes{
        {lay_sizes[1].first, lay_sizes[1].second}};
    for (unsigned int i = 2u; i < n_layers + 1u; ++i) {
        vol_sizes.emplace_back(lay_sizes[i].first, lay_sizes[i - 1u].second);
        vol_sizes.emplace_back(lay_sizes[i].first, lay_sizes[i].second);
    }

    brl_module_factory m_factory{cfg};
    bool is_gap = true;
    for (unsigned int i = 0u; i < 2u * n_layers - 1u; ++i) {
        unsigned int j = (i + 2u) / 2u;
        // Every second layer is a gap volume
        is_gap = !is_gap;
        if (is_gap) {
            detail::detector_helper<transform3_t>().create_cyl_volume(
                det, resource, ctx, vol_sizes[i].first, vol_sizes[i].second,
                -brl_half_z, brl_half_z, volume_links_vec[i]);

            dindex vol_idx = det.volumes().back().index();
            names[vol_idx + 1u] = "gap_" + std::to_string(vol_idx);
        } else {
            m_factory.cfg.m_binning = m_binning[j];
            m_factory.cfg.layer_r = lay_positions[j];
            detail::detector_helper<transform3_t>().create_cyl_volume(
                det, resource, ctx, vol_sizes[i].first, vol_sizes[i].second,
                -brl_half_z, brl_half_z, volume_links_vec[i]);

            dindex vol_idx = det.volumes().back().index();
            names[vol_idx + 1u] = "barrel_" + std::to_string(vol_idx);

            add_cylinder_grid(ctx, resource, det.volumes().back(), det, cfg,
                              m_factory);
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
    covfie::field<toy_metadata<>::bfield_backend_t> &&bfield,
    unsigned int n_brl_layers = 4u, unsigned int n_edc_layers = 3u) {

    // detector type
    using detector_t = detector<toy_metadata<>, covfie::field, container_t>;

    // Detector and volume names
    typename detector_t::name_map name_map = {{0u, "toy_detector"}};

    /// Leaving world
    using nav_link_t = typename detector_t::surface_type::navigation_link;
    constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};

    //
    // barrel
    //
    constexpr scalar brl_half_z{500.f};
    const std::vector<scalar> brl_positions = {19.f, 32.f, 72.f, 116.f, 172.f};
    const std::vector<std::pair<scalar, scalar>> brl_lay_sizes = {
        {0.f, 27.f},
        {27.f, 38.f},
        {64.f, 80.f},
        {108.f, 124.f},
        {164.f, 180.f}};
    const std::vector<std::pair<int, int>> brl_binning = {
        {0.f, 0.f}, {16.f, 14.f}, {32.f, 14.f}, {52.f, 14.f}, {78.f, 14.f}};
    // module parameters
    struct brl_m_config {
        scalar m_half_x{8.4f};
        scalar m_half_y{36.f};
        scalar m_tilt_phi{0.14f};  // 0.145;
        scalar layer_r{32.f};
        scalar m_radial_stagger{0.5f};  // 2.;
        scalar m_long_overlap{2.f};     // 5.;
        std::pair<unsigned int, unsigned int> m_binning = {16u, 14u};
        material<scalar> mat = silicon_tml<scalar>();
        scalar thickness{0.15f * unit<scalar>::mm};
    };

    //
    // endcaps
    //
    const std::vector<scalar> edc_positions = {600.f,  700.f,  820.f, 960.f,
                                               1100.f, 1300.f, 1500.f};
    const std::vector<std::pair<scalar, scalar>> edc_lay_sizes = {
        {595.f, 605.f},   {695.f, 705.f},   {815.f, 825.f},  {955.f, 965.f},
        {1095.f, 1105.f}, {1295.f, 1305.f}, {1495.f, 1505.f}};
    // module params
    struct edc_m_config {
        int side{1};
        scalar inner_r{27.f};
        scalar outer_r{180.f};
        scalar edc_position{600.f};
        scalar ring_stagger{1.f};
        // Parameters for both rings of modules
        std::vector<scalar> m_phi_stagger = {4.f, 4.f};
        std::vector<scalar> m_phi_sub_stagger = {0.5f, 0.5f};
        std::vector<unsigned int> disc_binning = {40u, 68u};
        std::vector<scalar> m_half_y = {36.f, 36.f};
        // std::vector<scalar> m_half_x_min_y = {8.4f, 8.4f};
        // std::vector<scalar> m_half_x_max_y = {10.1f, 10.1f};
        std::vector<scalar> m_half_x_min_y = {8.4f, 8.4f};
        std::vector<scalar> m_half_x_max_y = {12.4f, 12.4f};
        std::vector<scalar> m_tilt = {0.f, 0.f};
        material<scalar> mat = silicon_tml<scalar>();
        scalar thickness{0.15f * unit<scalar>::mm};
    };

    // Fills volume with barrel layer
    struct brl_module_factory {
        brl_m_config cfg;

        void operator()(const typename detector_t::geometry_context &ctx,
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

        void operator()(const typename detector_t::geometry_context &ctx,
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
    if (n_brl_layers > brl_positions.size() - 1u) {
        throw std::invalid_argument(
            "ERROR: Too many barrel layers requested (max " +
            std::to_string(brl_positions.size() - 1u) + ")!");
    }
    // the radius of the endcaps and  the barrel section need to match
    if (n_edc_layers > 0 and
        std::fabs(brl_lay_sizes[n_brl_layers].second - edc_config.outer_r) >
            std::numeric_limits<scalar>::epsilon()) {
        throw std::invalid_argument(
            "ERROR: Barrel and endcap radii do not match!");
    }

    // beampipe
    const dindex beampipe_idx{0u};
    add_beampipe(det, resource, ctx0, name_map, n_edc_layers, n_brl_layers,
                 edc_lay_sizes, brl_lay_sizes[0], brl_positions[0], brl_half_z,
                 edc_config.inner_r);

    if (n_edc_layers > 0u) {
        edc_config.side = -1;
        // negative endcap layers
        add_endcap_detector<edc_module_factory>(
            det, resource, ctx0, name_map, n_edc_layers, beampipe_idx,
            edc_lay_sizes, edc_positions, edc_config);

        // gap volume that connects barrel and neg. endcap
        dindex prev_vol_idx = det.volumes().back().index();
        prev_vol_idx = prev_vol_idx == 0u ? leaving_world : prev_vol_idx;
        dindex next_vol_idx = n_brl_layers == 0u
                                  ? leaving_world
                                  : det.volumes().back().index() + 2u;

        add_endcap_barrel_connection(
            det, resource, ctx0, name_map, edc_config.side, n_brl_layers,
            beampipe_idx, brl_lay_sizes, edc_config.inner_r, edc_config.outer_r,
            edc_lay_sizes[0].first, brl_half_z, next_vol_idx, prev_vol_idx);
    }
    if (n_brl_layers > 0u) {
        // barrel
        add_barrel_detector<brl_module_factory>(
            det, resource, ctx0, name_map, n_brl_layers, beampipe_idx,
            brl_half_z, brl_lay_sizes, brl_positions, brl_binning, brl_config);
    }
    if (n_edc_layers > 0u) {
        // gap layer that connects barrel and pos. endcap
        edc_config.side = 1.;
        // innermost barrel layer volume id
        dindex prev_vol_idx =
            n_brl_layers == 0 ? leaving_world : 2u * n_edc_layers + 1u;
        dindex next_vol_idx = prev_vol_idx == 1u
                                  ? leaving_world
                                  : det.volumes().back().index() + 2u;

        add_endcap_barrel_connection(
            det, resource, ctx0, name_map, edc_config.side, n_brl_layers,
            beampipe_idx, brl_lay_sizes, edc_config.inner_r, edc_config.outer_r,
            brl_half_z, edc_lay_sizes[0].first, prev_vol_idx, next_vol_idx);

        // positive endcap layers
        add_endcap_detector<edc_module_factory>(
            det, resource, ctx0, name_map, n_edc_layers, beampipe_idx,
            edc_lay_sizes, edc_positions, edc_config);
    }

    // Add volume grid
    // TODO: Fill it

    // Dimensions of the volume grid: minr, min phi, minz, maxr, maxphi, maxz
    // TODO: Adapt to number of layers
    mask<cylinder3D> vgrid_dims{0u,      0.f,   -constant<scalar>::pi,
                                -2000.f, 180.f, constant<scalar>::pi,
                                2000.f};
    std::array<std::size_t, 3> n_vgrid_bins{1u, 1u, 1u};

    grid_factory_type<typename detector_t::volume_finder> vgrid_factory{};
    auto vgrid =
        vgrid_factory.template new_grid<n_axis::open<n_axis::label::e_r>,
                                        n_axis::circular<n_axis::label::e_phi>,
                                        n_axis::open<n_axis::label::e_z>,
                                        n_axis::irregular<>, n_axis::regular<>,
                                        n_axis::irregular<>>(vgrid_dims,
                                                             n_vgrid_bins);
    det.set_volume_finder(std::move(vgrid));

    return std::make_pair(std::move(det), std::move(name_map));
}

/** Wrapper for create_toy_geometry with constant zero bfield.
 */
template <typename container_t = host_container_types>
auto create_toy_geometry(vecmem::memory_resource &resource,
                         unsigned int n_brl_layers = 4u,
                         unsigned int n_edc_layers = 3u) {
    return create_toy_geometry<container_t>(
        resource,
        covfie::field<toy_metadata<>::bfield_backend_t>{
            toy_metadata<>::bfield_backend_t::configuration_t{0.f, 0.f, 0.f}},
        n_brl_layers, n_edc_layers);
}

}  // namespace detray
