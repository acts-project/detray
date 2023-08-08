/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include <iostream>

#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/predefined_materials.hpp"

namespace detray::detail {

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
              typename mask_container_t, typename material_container_t,
              typename transform_container_t, typename volume_links>
    inline void add_cylinder_surface(
        const dindex volume_idx, context_t &ctx, surface_container_t &surfaces,
        mask_container_t &masks, material_container_t &materials,
        transform_container_t &transforms, const scalar_t r,
        const scalar_t lower_z, const scalar_t upper_z,
        const volume_links volume_link, const material<scalar_t> &mat,
        const scalar_t thickness) {
        using surface_type = typename surface_container_t::value_type;
        using mask_link_type = typename surface_type::mask_link;
        using material_id = typename surface_type::material_id;
        using material_link_type = typename surface_type::material_link;

        constexpr auto slab_id = material_id::e_slab;

        const scalar_t min_z{std::min(lower_z, upper_z)};
        const scalar_t max_z{std::max(lower_z, upper_z)};

        // translation
        point3 tsl{0.f, 0.f, 0.f};

        // add transform and masks
        transforms.emplace_back(ctx, tsl);
        masks.template emplace_back<cyl_id>(empty_context{}, volume_link, r,
                                            min_z, max_z);

        // Add material slab
        materials.template emplace_back<slab_id>(empty_context{}, mat,
                                                 thickness);

        // add surface
        mask_link_type mask_link{cyl_id, masks.template size<cyl_id>() - 1u};
        material_link_type material_link{
            slab_id, materials.template size<slab_id>() - 1u};
        const surface_id sf_id = (volume_link != volume_idx)
                                     ? surface_id::e_portal
                                     : surface_id::e_passive;

        surfaces.emplace_back(transforms.size(ctx) - 1u, mask_link,
                              material_link, volume_idx, dindex_invalid, sf_id);
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
        transform_container_t &transforms, const scalar_t inner_r,
        const scalar_t outer_r, const scalar_t z,
        const volume_links volume_link, const material<scalar_t> &mat,
        const scalar_t thickness) {
        using surface_type = typename surface_container_t::value_type;
        using mask_id = typename surface_type::mask_id;
        using mask_link_type = typename surface_type::mask_link;
        using material_id = typename surface_type::material_id;
        using material_link_type = typename surface_type::material_link;

        constexpr auto disc_id = mask_id::e_portal_ring2;
        constexpr auto slab_id = material_id::e_slab;

        const scalar_t min_r{std::min(inner_r, outer_r)};
        const scalar_t max_r{std::max(inner_r, outer_r)};

        // translation
        point3 tsl{0.f, 0.f, z};

        // add transform and mask
        transforms.emplace_back(ctx, tsl);
        masks.template emplace_back<disc_id>(empty_context{}, volume_link,
                                             min_r, max_r);

        // Add material slab
        materials.template emplace_back<slab_id>(empty_context{}, mat,
                                                 thickness);

        // add surface
        mask_link_type mask_link{disc_id, masks.template size<disc_id>() - 1u};
        material_link_type material_link{
            slab_id, materials.template size<slab_id>() - 1u};
        const surface_id sf_id = (volume_link != volume_idx)
                                     ? surface_id::e_portal
                                     : surface_id::e_sensitive;
        surfaces.emplace_back(transforms.size(ctx) - 1u, mask_link,
                              material_link, volume_idx, dindex_invalid, sf_id);
    }

    /** Function that adds a generic cylinder volume, using a factory for
     * contained module surfaces.
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
    template <typename detector_t>
    void create_cyl_volume(detector_t &det, vecmem::memory_resource &resource,
                           typename detector_t::geometry_context &ctx,
                           const scalar_t lay_inner_r,
                           const scalar_t lay_outer_r, const scalar_t lay_neg_z,
                           const scalar_t lay_pos_z,
                           const std::vector<dindex> &volume_links) {
        // volume bounds
        const scalar_t inner_r{std::min(lay_inner_r, lay_outer_r)};
        const scalar_t outer_r{std::max(lay_inner_r, lay_outer_r)};
        const scalar_t lower_z{std::min(lay_neg_z, lay_pos_z)};
        const scalar_t upper_z{std::max(lay_neg_z, lay_pos_z)};

        // Add module surfaces to volume
        typename detector_t::surface_container_t surfaces(&resource);
        typename detector_t::mask_container masks(resource);
        typename detector_t::material_container materials(resource);
        typename detector_t::transform_container transforms(resource);

        auto &cyl_volume = det.new_volume(
            volume_id::e_cylinder, {detector_t::sf_finders::id::e_default, 0u});

        // volume placement
        cyl_volume.set_transform(det.transform_store().size());
        // translation of the cylinder
        point3 t{0.f, 0.f, 0.5f * (upper_z + lower_z)};
        det.transform_store().emplace_back(ctx, t);

        // negative and positive, inner and outer portal surface
        constexpr auto cyl_id = detector_t::masks::id::e_portal_cylinder2;

        // If inner radius is 0, skip add the innder cylinder
        if (inner_r > 0.f) {
            add_cylinder_surface<cyl_id>(
                cyl_volume.index(), ctx, surfaces, masks, materials, transforms,
                inner_r, lower_z, upper_z, volume_links[0], vacuum<scalar_t>(),
                0.f * unit<scalar_t>::mm);
        }
        add_cylinder_surface<cyl_id>(
            cyl_volume.index(), ctx, surfaces, masks, materials, transforms,
            outer_r, lower_z, upper_z, volume_links[1], vacuum<scalar_t>(),
            0.f * unit<scalar_t>::mm);
        add_disc_surface(cyl_volume.index(), ctx, surfaces, masks, materials,
                         transforms, inner_r, outer_r, lower_z, volume_links[2],
                         vacuum<scalar_t>(), 0.f * unit<scalar_t>::mm);
        add_disc_surface(cyl_volume.index(), ctx, surfaces, masks, materials,
                         transforms, inner_r, outer_r, upper_z, volume_links[3],
                         vacuum<scalar_t>(), 0.f * unit<scalar_t>::mm);

        det.add_objects_per_volume(ctx, cyl_volume, surfaces, masks, transforms,
                                   materials);
    }
};

}  // namespace detray::detail