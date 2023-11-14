/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/materials/material_map.hpp"
#include "detray/materials/material_slab.hpp"

namespace detray::detail {

/// Generate material along z bins for a cylinder material grid
inline std::vector<material_slab<scalar>> generate_cyl_mat(
    const std::array<scalar, 2u> &bounds, const std::size_t nbins,
    material<scalar> mat, const scalar t) {
    std::vector<material_slab<scalar>> ts;
    ts.reserve(nbins);

    scalar z{bounds[0]};
    const scalar z_step{(bounds[1] - bounds[0]) /
                        static_cast<scalar>(nbins - 1u)};
    for (std::size_t n = 0u; n < nbins; ++n) {
        ts.emplace_back(mat, static_cast<scalar>(0.00001f * z * z) + t);
        z += z_step;
    }

    return ts;
}

/// Generate material along r bins for a disc material grid
inline std::vector<material_slab<scalar>> generate_disc_mat(
    const std::array<scalar, 2u> &bounds, const std::size_t nbins,
    material<scalar> mat, const scalar t) {
    std::vector<material_slab<scalar>> ts;
    ts.reserve(nbins);

    scalar r{bounds[0]};
    const scalar r_step{(bounds[1] - bounds[0]) /
                        static_cast<scalar>(nbins - 1u)};
    for (std::size_t n = 0u; n < nbins; ++n) {
        ts.emplace_back(mat, static_cast<scalar>(0.01f * r) + t);
        r += r_step;
    }

    return ts;
}

template <typename algebra_t>
struct detector_helper {

    using scalar_t = typename algebra_t::scalar_type;
    using point3 = typename algebra_t::point3;

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
              typename volume_links>
    inline auto add_cylinder_surface(const dindex volume_idx, context_t &ctx,
                                     surface_container_t &surfaces,
                                     mask_container_t &masks,
                                     transform_container_t &transforms,
                                     const scalar_t r, const scalar_t lower_z,
                                     const scalar_t upper_z,
                                     const volume_links volume_link) const {
        using surface_type = typename surface_container_t::value_type;
        using mask_link_type = typename surface_type::mask_link;
        using material_id = typename surface_type::material_id;
        using material_link_type = typename surface_type::material_link;

        const scalar_t min_z{math::min(lower_z, upper_z)};
        const scalar_t max_z{math::max(lower_z, upper_z)};

        // translation
        point3 tsl{0.f, 0.f, 0.f};

        // add transform and masks
        transforms.emplace_back(ctx, tsl);
        auto &mask_ref = masks.template emplace_back<cyl_id>(
            empty_context{}, volume_link, r, min_z, max_z);

        // add surface
        mask_link_type mask_link{cyl_id, masks.template size<cyl_id>() - 1u};
        material_link_type material_link{material_id::e_none, dindex_invalid};
        const surface_id sf_id = (volume_link != volume_idx)
                                     ? surface_id::e_portal
                                     : surface_id::e_passive;

        auto &sf_ref =
            surfaces.emplace_back(transforms.size(ctx) - 1u, mask_link,
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
              typename volume_links>
    inline auto add_disc_surface(const dindex volume_idx, context_t &ctx,
                                 surface_container_t &surfaces,
                                 mask_container_t &masks,
                                 transform_container_t &transforms,
                                 const scalar_t inner_r, const scalar_t outer_r,
                                 const scalar_t z,
                                 const volume_links volume_link) const {
        using surface_type = typename surface_container_t::value_type;
        using mask_id = typename surface_type::mask_id;
        using mask_link_type = typename surface_type::mask_link;
        using material_id = typename surface_type::material_id;
        using material_link_type = typename surface_type::material_link;

        constexpr auto disc_id = mask_id::e_portal_ring2;

        const scalar_t min_r{math::min(inner_r, outer_r)};
        const scalar_t max_r{math::max(inner_r, outer_r)};

        // translation
        point3 tsl{0.f, 0.f, z};

        // add transform and mask
        transforms.emplace_back(ctx, tsl);
        auto &mask_ref = masks.template emplace_back<disc_id>(
            empty_context{}, volume_link, min_r, max_r);

        // add surface
        mask_link_type mask_link{disc_id, masks.template size<disc_id>() - 1u};
        material_link_type material_link{material_id::e_none, dindex_invalid};
        const surface_id sf_id = (volume_link != volume_idx)
                                     ? surface_id::e_portal
                                     : surface_id::e_sensitive;
        auto &sf_ref =
            surfaces.emplace_back(transforms.size(ctx) - 1u, mask_link,
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
    void create_cyl_volume(const config_t &cfg, detector_t &det,
                           vecmem::memory_resource &resource,
                           typename detector_t::geometry_context &ctx,
                           const scalar_t lay_inner_r,
                           const scalar_t lay_outer_r, const scalar_t lay_neg_z,
                           const scalar_t lay_pos_z,
                           const std::vector<dindex> &volume_links) const {
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

        auto &cyl_volume = det.new_volume(
            volume_id::e_cylinder, {detector_t::accel::id::e_default, 0u});

        // volume placement
        cyl_volume.set_transform(det.transform_store().size());
        // translation of the cylinder
        point3 t{0.f, 0.f, 0.5f * (upper_z + lower_z)};
        det.transform_store().emplace_back(ctx, t);

        // negative and positive, inner and outer portal surface
        constexpr auto cyl_id = detector_t::masks::id::e_portal_cylinder2;

        auto &material_coll =
            cfg.use_material_maps() ? det.material_store() : materials;

        // If inner radius is 0, skip adding the inner cylinder
        if (inner_r > 0.f) {
            auto [inner_cyl, inner_cyl_mask] = add_cylinder_surface<cyl_id>(
                cyl_volume.index(), ctx, surfaces, masks, transforms, inner_r,
                lower_z, upper_z, volume_links[0]);
            create_material(cfg, inner_cyl, inner_cyl_mask, material_coll);
        }

        auto [outer_cyl, outer_cyl_mask] = add_cylinder_surface<cyl_id>(
            cyl_volume.index(), ctx, surfaces, masks, transforms, outer_r,
            lower_z, upper_z, volume_links[1]);
        create_material(cfg, outer_cyl, outer_cyl_mask, material_coll);

        auto [neg_disc, neg_disc_mask] = add_disc_surface(
            cyl_volume.index(), ctx, surfaces, masks, transforms, inner_r,
            outer_r, lower_z, volume_links[2]);
        create_material(cfg, neg_disc, neg_disc_mask, material_coll);

        auto [pos_disc, pos_disc_mask] = add_disc_surface(
            cyl_volume.index(), ctx, surfaces, masks, transforms, inner_r,
            outer_r, upper_z, volume_links[3]);
        create_material(cfg, pos_disc, pos_disc_mask, material_coll);

        det.add_objects_per_volume(ctx, cyl_volume, surfaces, masks, transforms,
                                   materials);

        // Surface index offset in the global detector container
        auto sf_offset{static_cast<dindex>(det.surfaces().size())};

        // Identify the portal index range in the global container
        dindex_range pt_range{sf_offset - (inner_r > 0.f ? 4u : 3u), sf_offset};

        det.volumes().back().template update_sf_link<surface_id::e_portal>(
            pt_range);
    }

    template <typename config_t, typename surface_desc_t, typename mask_t,
              typename material_container_t>
    inline void create_material(const config_t &cfg, surface_desc_t &sf,
                                const mask_t &sf_mask,
                                material_container_t &materials) const {
        using material_id = typename surface_desc_t::material_id;
        using material_link_type = typename surface_desc_t::material_link;

        if (cfg.use_material_maps()) {
            material_grid_factory<scalar_t> mat_map_factory{};

            constexpr auto is_disc_map{
                std::is_same_v<typename mask_t::shape, ring2D<>>};

            // Scale material thickness either over r- or z-bins
            const std::size_t bins =
                is_disc_map ? cfg.disc_map_bins()[0] : cfg.cyl_map_bins()[1];

            const auto &bounds = sf_mask.values();
            const auto mat{
                is_disc_map
                    ? cfg.edc_mat_generator()({bounds[ring2D<>::e_inner_r],
                                               bounds[ring2D<>::e_outer_r]},
                                              bins, cfg.mapped_material(),
                                              cfg.thickness())
                    : cfg.barrel_mat_generator()(
                          {bounds[cylinder2D<>::e_n_half_z],
                           bounds[cylinder2D<>::e_p_half_z]},
                          bins, cfg.mapped_material(), cfg.thickness())};

            auto material_map = mat_map_factory.template new_grid<>(
                sf_mask,
                is_disc_map ? cfg.disc_map_bins() : cfg.cyl_map_bins());

            const auto n_bins_per_axis{material_map.axes().nbins_per_axis()};
            for (dindex bin0 = 0u; bin0 < n_bins_per_axis[0]; ++bin0) {
                for (dindex bin1 = 0u; bin1 < n_bins_per_axis[1]; ++bin1) {
                    material_map.template populate<replace<>>(
                        axis::multi_bin<2>{bin0, bin1},
                        mat[is_disc_map ? bin0 : bin1]);
                }
            }

            constexpr auto map_id{is_disc_map ? material_id::e_disc2_map
                                              : material_id::e_cylinder2_map};

            materials.template push_back<map_id>(material_map);

            sf.material() = material_link_type{
                map_id, materials.template size<map_id>() - 1u};
        } else {
            materials.template emplace_back<material_id::e_slab>(
                {}, cfg.mapped_material(), cfg.thickness());

            sf.material() = material_link_type{
                material_id::e_slab,
                materials.template size<material_id::e_slab>() - 1u};
        }
    }
};

}  // namespace detray::detail
