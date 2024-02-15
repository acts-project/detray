/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/indexing.hpp"
#include "detray/materials/predefined_materials.hpp"

namespace detray {

/// Adds a few surfaces to the detector.
template <typename detector_t>
void prefill_detector(detector_t& d,
                      typename detector_t::geometry_context ctx) {
    using scalar_t = typename detector_t::scalar_type;
    using point3 = typename detector_t::point3;
    using mask_id = typename detector_t::masks::id;
    using material_id = typename detector_t::materials::id;
    using surface_t = typename detector_t::surface_type;
    using mask_link_t = typename surface_t::mask_link;
    using material_link_t = typename surface_t::material_link;

    detray::empty_context empty_ctx{};
    vecmem::memory_resource* host_mr = d.resource();
    typename detector_t::transform_container trfs(*host_mr);
    typename detector_t::surface_container surfaces{};
    typename detector_t::mask_container masks(*host_mr);
    typename detector_t::material_container materials(*host_mr);

    /// Surface 0
    point3 t0{0.f, 0.f, 0.f};
    trfs.emplace_back(ctx, t0);
    masks.template emplace_back<mask_id::e_rectangle2>(empty_ctx, 0u, -3.f,
                                                       3.f);
    materials.template emplace_back<material_id::e_slab>(
        empty_ctx, detray::gold<scalar_t>(), 3.f);
    mask_link_t mask_link{mask_id::e_rectangle2,
                          masks.template size<mask_id::e_rectangle2>() - 1};
    material_link_t material_link{
        material_id::e_slab,
        materials.template size<material_id::e_slab>() - 1};
    surfaces.emplace_back(trfs.size(ctx) - 1, mask_link, material_link, 0u,
                          detray::surface_id::e_sensitive);
    surfaces.back().set_index(
        static_cast<detray::dindex>(surfaces.size() - 1u));

    /// Surface 1
    point3 t1{1.f, 0.f, 0.f};
    trfs.emplace_back(ctx, t1);
    masks.template emplace_back<mask_id::e_annulus2>(empty_ctx, 0u, 1.f, 2.f,
                                                     3.f, 4.f, 5.f, 6.f, 7.f);
    materials.template emplace_back<material_id::e_slab>(
        empty_ctx, detray::tungsten<scalar_t>(), 12.f);

    mask_link = {mask_id::e_annulus2,
                 masks.template size<mask_id::e_annulus2>() - 1};
    material_link = {material_id::e_slab,
                     materials.template size<material_id::e_slab>() - 1};
    surfaces.emplace_back(trfs.size(ctx) - 1, mask_link, material_link, 0u,
                          detray::surface_id::e_sensitive);
    surfaces.back().set_index(
        static_cast<detray::dindex>(surfaces.size() - 1u));

    /// Surface 2
    point3 t2{2.f, 0.f, 0.f};
    trfs.emplace_back(ctx, t2);
    masks.template emplace_back<mask_id::e_trapezoid2>(empty_ctx, 0u, 1.f, 2.f,
                                                       3.f);
    materials.template emplace_back<material_id::e_rod>(
        empty_ctx, detray::aluminium<scalar_t>(), 4.f);

    mask_link = {mask_id::e_trapezoid2,
                 masks.template size<mask_id::e_trapezoid2>() - 1};
    material_link = {material_id::e_rod,
                     materials.template size<material_id::e_rod>() - 1};
    surfaces.emplace_back(trfs.size(ctx) - 1, mask_link, material_link, 0u,
                          detray::surface_id::e_sensitive);
    surfaces.back().set_index(
        static_cast<detray::dindex>(surfaces.size() - 1u));

    // Add surfaces to lookup, so they can be easily fetched using a barcode
    for (const auto sf : surfaces) {
        d.surfaces().insert(sf);
    }

    // Add the new data
    d.new_volume(detray::volume_id::e_cylinder);
    d.append_portals(std::move(surfaces));
    d.append_transforms(std::move(trfs));
    d.append_masks(std::move(masks));
    d.append_materials(std::move(materials));
}

}  // namespace detray
