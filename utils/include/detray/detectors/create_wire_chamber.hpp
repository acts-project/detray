/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/grid_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/core/detector_metadata.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/materials/material_map.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/materials/mixture.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/utils/axis_rotation.hpp"
#include "detray/utils/unit_vectors.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

namespace detray {

struct wire_chamber_config {

    /// Number of layers
    unsigned int m_n_layers{10u};
    /// Half z of cylinder chamber
    scalar m_half_z{1000.f * unit<scalar>::mm};

    constexpr wire_chamber_config &n_layers(const unsigned int n) {
        m_n_layers = n;
        return *this;
    }
    constexpr wire_chamber_config &half_z(const scalar hz) {
        m_half_z = hz;
        return *this;
    }

    constexpr unsigned int n_layers() const { return m_n_layers; }
    constexpr scalar half_z() const { return m_half_z; }

};  // wire chamber config

namespace detail {

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
template <auto cyl_id, typename context_t, typename surface_container_t,
          typename mask_container_t, typename transform_container_t,
          typename volume_links, typename scalar_t>
inline auto add_cylinder_surface(const dindex volume_idx, context_t &ctx,
                                 surface_container_t &surfaces,
                                 mask_container_t &masks,
                                 transform_container_t &transforms,
                                 const scalar_t r, const scalar_t lower_z,
                                 const scalar_t upper_z,
                                 const volume_links volume_link) {
    using surface_type = typename surface_container_t::value_type;
    using mask_link_type = typename surface_type::mask_link;
    using material_id = typename surface_type::material_id;
    using material_link_type = typename surface_type::material_link;

    const scalar_t min_z{math::min(lower_z, upper_z)};
    const scalar_t max_z{math::max(lower_z, upper_z)};

    // translation
    typename transform_container_t::value_type::point3 tsl{0.f, 0.f, 0.f};

    // add transform and masks
    transforms.emplace_back(ctx, tsl);
    auto &mask_ref = masks.template emplace_back<cyl_id>(
        empty_context{}, volume_link, r, min_z, max_z);

    // add surface
    mask_link_type mask_link{cyl_id, masks.template size<cyl_id>() - 1u};
    material_link_type material_link{material_id::e_none, 0u};
    const surface_id sf_id = (volume_link != volume_idx)
                                 ? surface_id::e_portal
                                 : surface_id::e_passive;

    auto &sf_ref = surfaces.emplace_back(transforms.size(ctx) - 1u, mask_link,
                                         material_link, volume_idx, sf_id);

    return std::tie(sf_ref, mask_ref);
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
 * @param inner_r lower radius of the disc
 * @param outer_r upper radius of the disc
 * @param z z position of the disc
 * @param volume_link link to next volume for the masks
 */
template <typename context_t, typename surface_container_t,
          typename mask_container_t, typename transform_container_t,
          typename volume_links, typename scalar_t>
inline auto add_disc_surface(const dindex volume_idx, context_t &ctx,
                             surface_container_t &surfaces,
                             mask_container_t &masks,
                             transform_container_t &transforms,
                             const scalar_t inner_r, const scalar_t outer_r,
                             const scalar_t z, const volume_links volume_link) {
    using surface_type = typename surface_container_t::value_type;
    using mask_id = typename surface_type::mask_id;
    using mask_link_type = typename surface_type::mask_link;
    using material_id = typename surface_type::material_id;
    using material_link_type = typename surface_type::material_link;

    constexpr auto disc_id = mask_id::e_portal_ring2;

    const scalar_t min_r{math::min(inner_r, outer_r)};
    const scalar_t max_r{math::max(inner_r, outer_r)};

    // translation
    typename transform_container_t::value_type::point3 tsl{0.f, 0.f, z};

    // add transform and mask
    transforms.emplace_back(ctx, tsl);
    auto &mask_ref = masks.template emplace_back<disc_id>(
        empty_context{}, volume_link, min_r, max_r);

    // add surface
    mask_link_type mask_link{disc_id, masks.template size<disc_id>() - 1u};
    material_link_type material_link{material_id::e_none, 0u};
    const surface_id sf_id = (volume_link != volume_idx)
                                 ? surface_id::e_portal
                                 : surface_id::e_sensitive;
    auto &sf_ref = surfaces.emplace_back(transforms.size(ctx) - 1u, mask_link,
                                         material_link, volume_idx, sf_id);

    return std::tie(sf_ref, mask_ref);
}

/** Function that adds a generic cylinder volume, using a factory for
 * contained module surfaces.
 *
 * @tparam detector_t the detector type
 * @tparam factory_t type of module factory. Must be callable on containers.
 *
 * @param cfg detector building configuration
 * @param det detector the volume should be added to
 * @param resource vecmem memory resource for the temporary containers
 * @param ctx geometry context
 * @param lay_inner_r inner radius of volume
 * @param lay_outer_r outer radius of volume
 * @param lay_neg_r lower extend of volume
 * @param lay_pos_r upper extend of volume
 * @param volume_links volume links for the portals of the volume
 */
template <typename config_t, typename detector_t>
void create_cyl_volume(const config_t & /*cfg*/, detector_t &det,
                       vecmem::memory_resource &resource,
                       typename detector_t::geometry_context &ctx,
                       const typename detector_t::scalar_type lay_inner_r,
                       const typename detector_t::scalar_type lay_outer_r,
                       const typename detector_t::scalar_type lay_neg_z,
                       const typename detector_t::scalar_type lay_pos_z,
                       const std::vector<dindex> &volume_links) {

    using scalar_t = typename detector_t::scalar_type;
    using point3_t = typename detector_t::point3_type;

    // volume bounds
    const scalar_t inner_r{math::min(lay_inner_r, lay_outer_r)};
    const scalar_t outer_r{math::max(lay_inner_r, lay_outer_r)};
    const scalar_t lower_z{math::min(lay_neg_z, lay_pos_z)};
    const scalar_t upper_z{math::max(lay_neg_z, lay_pos_z)};

    // Add module surfaces to volume
    typename detector_t::surface_container surfaces(&resource);
    typename detector_t::mask_container masks(resource);
    typename detector_t::material_container materials(resource);
    typename detector_t::transform_container transforms(resource);

    auto &cyl_volume = det.new_volume(volume_id::e_cylinder,
                                      {detector_t::accel::id::e_default, 0u});
    cyl_volume.set_material(detector_t::volume_type::material_id::e_none, 0u);

    // volume placement
    cyl_volume.set_transform(det.transform_store().size());
    // translation of the cylinder
    point3_t t{0.f, 0.f, 0.5f * (upper_z + lower_z)};
    det.transform_store().emplace_back(ctx, t);

    // negative and positive, inner and outer portal surface
    constexpr auto cyl_id = detector_t::masks::id::e_portal_cylinder2;

    // If inner radius is 0, skip adding the inner cylinder
    if (inner_r > 0.f) {
        add_cylinder_surface<cyl_id>(cyl_volume.index(), ctx, surfaces, masks,
                                     transforms, inner_r, lower_z, upper_z,
                                     volume_links[0]);
    }

    add_cylinder_surface<cyl_id>(cyl_volume.index(), ctx, surfaces, masks,
                                 transforms, outer_r, lower_z, upper_z,
                                 volume_links[1]);

    add_disc_surface(cyl_volume.index(), ctx, surfaces, masks, transforms,
                     inner_r, outer_r, lower_z, volume_links[2]);

    add_disc_surface(cyl_volume.index(), ctx, surfaces, masks, transforms,
                     inner_r, outer_r, upper_z, volume_links[3]);

    det.add_objects_per_volume(ctx, cyl_volume, surfaces, masks, transforms,
                               materials);

    // Surface index offset in the global detector container
    auto sf_offset{static_cast<dindex>(det.surfaces().size())};

    // Identify the portal index range in the global container
    dindex_range pt_range{sf_offset - (inner_r > 0.f ? 4u : 3u), sf_offset};

    det.volumes().back().template update_sf_link<surface_id::e_portal>(
        pt_range);
}

}  // namespace detail

inline auto create_wire_chamber(vecmem::memory_resource &resource,
                                const wire_chamber_config &cfg) {

    // Detector type
    using detector_t = detector<default_metadata, host_container_types>;

    using algebra_t = detector_t::algebra_type;
    using scalar_t = typename detector_t::scalar_type;
    using point3_t = detector_t::point3_type;
    using vector3_t = detector_t::vector3_type;

    using nav_link_t = typename detector_t::surface_type::navigation_link;
    using mask_id = typename detector_t::surface_type::mask_id;
    using material_id = typename detector_t::surface_type::material_id;
    using mask_link_type = typename detector_t::surface_type::mask_link;
    using material_link_type = typename detector_t::surface_type::material_link;
    constexpr auto wire_id = mask_id::e_drift_cell;
    constexpr auto rod_id = material_id::e_rod;
    constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};

    // Detector configurations
    const scalar_t cyl_half_z{cfg.half_z()};
    constexpr scalar_t inner_cyl_rad{500.f * unit<scalar_t>::mm};
    constexpr scalar_t cell_size = 10.f * unit<scalar_t>::mm;
    constexpr scalar_t stereo_angle = 50.f * unit<scalar_t>::mrad;
    const material<scalar_t> wire_mat = tungsten<scalar_t>();
    constexpr scalar_t wire_rad = 15.f * unit<scalar_t>::um;

    // Create detector
    detector_t det(resource);

    // Detector and volume names
    typename detector_t::name_map name_map = {{0u, "wire_chamber"}};

    // geometry context object
    typename detector_t::geometry_context ctx0{};

    // Beam collision volume
    detail::create_cyl_volume(
        cfg, det, resource, ctx0, 0.f, inner_cyl_rad, -cyl_half_z, cyl_half_z,
        {leaving_world, 1u, leaving_world, leaving_world});

    name_map[1u] = "beam_vol_0";

    // Layer volumes
    const unsigned int n_layers{cfg.n_layers()};
    for (unsigned int i_lay = 0; i_lay < n_layers; i_lay++) {

        // Create a volume
        const scalar_t inner_layer_rad =
            inner_cyl_rad + static_cast<scalar_t>(i_lay) * cell_size * 2.f;
        const scalar_t outer_layer_rad =
            inner_cyl_rad + static_cast<scalar_t>(i_lay + 1) * cell_size * 2.f;

        if (i_lay < n_layers - 1) {
            detail::create_cyl_volume(
                cfg, det, resource, ctx0, inner_layer_rad, outer_layer_rad,
                -cyl_half_z, cyl_half_z,
                {i_lay, i_lay + 2, leaving_world, leaving_world});
        } else {
            detail::create_cyl_volume(
                cfg, det, resource, ctx0, inner_layer_rad, outer_layer_rad,
                -cyl_half_z, cyl_half_z,
                {i_lay, leaving_world, leaving_world, leaving_world});
        }

        // Current vol
        auto &vol_desc = det.volumes().back();

        // Layer configuration
        const scalar_t center_layer_rad = inner_layer_rad + cell_size;
        const scalar_t delta = 2 * cell_size / center_layer_rad;

        // Get volume ID
        auto volume_idx = vol_desc.index();
        name_map[volume_idx + 1u] = "layer_vol_" + std::to_string(volume_idx);

        auto mask_volume_link{static_cast<nav_link_t>(volume_idx)};

        // Containers per volume
        typename detector_t::surface_container surfaces(&resource);
        typename detector_t::mask_container masks(resource);
        typename detector_t::material_container materials(resource);
        typename detector_t::transform_container transforms(resource);

        // Wire center positions
        detray::dvector<point3_t> m_centers{};

        scalar_t theta{0.f};
        while (theta <= 2.f * constant<scalar_t>::pi) {

            const scalar_t x = center_layer_rad * math::cos(theta);
            const scalar_t y = center_layer_rad * math::sin(theta);
            const scalar_t z = 0.f;

            m_centers.push_back({x, y, z});
            theta += delta;
        }

        for (auto &m_center : m_centers) {

            // Surfaces with the linking into the local containers
            mask_link_type mask_link = {wire_id,
                                        masks.template size<wire_id>()};
            material_link_type material_link{rod_id,
                                             materials.template size<rod_id>()};
            const auto trf_index = transforms.size(ctx0);
            surfaces.emplace_back(trf_index, mask_link, material_link,
                                  volume_idx, surface_id::e_sensitive);

            // The wire bounds
            masks.template emplace_back<wire_id>(
                empty_context{}, mask_volume_link, cell_size, cyl_half_z);
            materials.template emplace_back<rod_id>(empty_context{}, wire_mat,
                                                    wire_rad);

            // Build the transform
            vector3_t z_axis{0.f, 0.f, 1.f};
            vector3_t r_axis = vector::normalize(m_center);
            const scalar sign = (i_lay % 2 == 0) ? 1 : -1;
            z_axis =
                axis_rotation<algebra_t>(r_axis, sign * stereo_angle)(z_axis);
            vector3_t x_axis =
                unit_vectors<vector3_t>().make_curvilinear_unit_u(z_axis);
            transforms.emplace_back(ctx0, m_center, z_axis, x_axis);
        }

        // Iterate the surfaces and update their links
        const auto trf_offset{det.transform_store().size(ctx0)};
        auto sf_offset{static_cast<dindex>(det.surfaces().size())};

        vol_desc.template update_sf_link<surface_id::e_sensitive>(
            sf_offset, surfaces.size());

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
            det.surfaces().insert(sf_desc);
        }

        // Add transforms, masks and material to detector
        det.append_masks(std::move(masks));
        det.append_transforms(std::move(transforms));
        det.append_materials(std::move(materials));

        //
        // Fill Grid
        //

        // Get relevant ids
        using geo_obj_ids = typename detector_t::geo_obj_ids;
        constexpr auto cyl_id = detector_t::masks::id::e_portal_cylinder2;
        constexpr auto grid_id = detector_t::accel::id::e_cylinder2_grid;

        using cyl_grid_t =
            typename detector_t::accelerator_container::template get_type<
                grid_id>;
        auto gbuilder = grid_builder<detector_t, cyl_grid_t>{};

        // The portal portals are at the end of the portal range by construction
        auto vol = detector_volume{det, vol_desc};
        auto portal_mask_idx = (vol.portals().end() - 4)->mask().index();
        const auto &inner_cyl_mask =
            det.mask_store().template get<cyl_id>().at(portal_mask_idx);

        portal_mask_idx = (vol.portals().end() - 3)->mask().index();
        const auto &outer_cyl_mask =
            det.mask_store().template get<cyl_id>().at(portal_mask_idx);

        // Correct cylinder radius so that the grid lies in the middle
        using cyl_mask_t = detail::remove_cvref_t<decltype(outer_cyl_mask)>;
        typename cyl_mask_t::mask_values mask_values{outer_cyl_mask.values()};
        mask_values[cylinder2D::e_r] =
            0.5f * (inner_cyl_mask.values()[cylinder2D::e_r] +
                    outer_cyl_mask.values()[cylinder2D::e_r]);
        const cyl_mask_t cyl_mask{mask_values, 0u};

        std::vector<std::pair<typename cyl_grid_t::loc_bin_index, dindex>>
            capacities{};
        auto bin_indexer2D = detray::views::cartesian_product{
            detray::views::iota{0u, 100u}, detray::views::iota{0u, 1u}};
        for (const auto &bin_idx : bin_indexer2D) {
            typename cyl_grid_t::loc_bin_index mbin{std::get<0>(bin_idx),
                                                    std::get<1>(bin_idx)};
            // @Todo: fine-tune capacity
            capacities.emplace_back(mbin, 3u);
        }

        // Add new grid to the detector
        gbuilder.init_grid(cyl_mask, {100u, 1u}, capacities);
        gbuilder.fill_grid(vol, det.surfaces(), det.transform_store(),
                           det.mask_store(), ctx0);

        det.accelerator_store().template push_back<grid_id>(gbuilder.get());
        vol_desc.template set_accel_link<geo_obj_ids::e_sensitive>(
            grid_id, det.accelerator_store().template size<grid_id>() - 1u);

        // Add volume grid
        // TODO: Fill it

        // Dimensions of the volume grid: minr, min phi, minz, maxr, maxphi,
        // maxz
        // TODO: Adapt to number of layers
        mask<cylinder3D> vgrid_dims{0u,     0.f,   -constant<scalar_t>::pi,
                                    -600.f, 180.f, constant<scalar_t>::pi,
                                    600.f};
        std::array<std::size_t, 3> n_vgrid_bins{1u, 1u, 1u};

        const scalar_t outer_radius{inner_cyl_rad +
                                    static_cast<scalar_t>(cfg.n_layers() + 1) *
                                        cell_size * 2.f};
        std::array<std::vector<scalar_t>, 3UL> bin_edges{
            std::vector<scalar_t>{0.f, outer_radius},
            std::vector<scalar_t>{-constant<scalar_t>::pi,
                                  constant<scalar_t>::pi},
            std::vector<scalar_t>{-cyl_half_z, cyl_half_z}};

        grid_factory_type<typename detector_t::volume_finder> vgrid_factory{};
        auto vgrid = vgrid_factory.template new_grid<
            axis::open<axis::label::e_r>, axis::circular<axis::label::e_phi>,
            axis::open<axis::label::e_z>, axis::irregular<>, axis::regular<>,
            axis::irregular<>>(vgrid_dims, n_vgrid_bins, {}, bin_edges);
        det.set_volume_finder(std::move(vgrid));
    }

    return std::make_pair(std::move(det), std::move(name_map));
}

}  // namespace detray
