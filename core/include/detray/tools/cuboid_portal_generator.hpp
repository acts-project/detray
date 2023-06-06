/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/masks/cuboid3D.hpp"
#include "detray/tools/bounding_volume.hpp"
#include "detray/tools/surface_factory_interface.hpp"

// System include(s)
#include <cassert>
#include <limits>

namespace detray {

/// @brief Generates a portal box around a volume that already contains surfaces
///
/// @tparam detector_t the type of detector the volume belongs to.
template <typename detector_t>
class cuboid_portal_generator final
    : public surface_factory_interface<detector_t> {

    using scalar_t = typename detector_t::scalar_type;
    using transform3_t = typename detector_t::transform3;

    /// A functor to construct global bounding boxes around masks
    struct bounding_box_creator {

        using aabb_t = axis_aligned_bounding_volume<cuboid3D<>>;

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline void operator()(
            const mask_group_t &mask_group, const index_t &index,
            const scalar_t envelope, const transform3_t &trf,
            std::vector<aabb_t> &boxes) const {
            // Local minimum bounding box
            aabb_t box{mask_group[index], boxes.size(), envelope};
            // Bounding box in global coordinates (might no longer be minimum)
            boxes.push_back(box.transform(trf));
        }
    };

    public:
    /// Use @param env as portal envelope
    DETRAY_HOST
    cuboid_portal_generator(const scalar_t env) : m_envelope{env} {}

    /// @returns the id for the surface type (portal surfaces)
    DETRAY_HOST
    auto surface_type() const -> surface_id override {
        return surface_id::e_portal;
    }

    /// @returns the number of rectangle portals this factory will produce
    DETRAY_HOST
    auto size() const -> dindex override { return 6u; }

    DETRAY_HOST
    void push_back(surface_data<detector_t> &&) override { /*Do nothing*/
    }
    DETRAY_HOST
    auto push_back(std::vector<surface_data<detector_t>> &&)
        -> void override { /*Do nothing*/
    }

    /// Create minimum aabbs around all surfaces that are passed and then
    /// construct world portals from the global bounding box.
    ///
    /// @param volume the volume the portals need to be added to.
    /// @param surfaces the surface collection to wrap and to add the portals to
    /// @param transforms the transforms of the surfaces.
    /// @param masks the masks of the surfaces.
    /// @param ctx the geometry context (not needed for portals).
    DETRAY_HOST
    auto operator()(typename detector_t::volume_type &volume,
                    typename detector_t::surface_container_t &surfaces,
                    typename detector_t::transform_container &transforms,
                    typename detector_t::mask_container &masks,
                    typename detector_t::geometry_context ctx = {}) const
        -> dindex_range override {

        using point3_t = typename detector_t::point3;
        using vector3_t = typename detector_t::vector3;

        using surface_t = typename detector_t::surface_type;
        using nav_link_t = typename surface_t::navigation_link;
        using mask_link_t = typename surface_t::mask_link;
        using material_link_t = typename surface_t::material_link;

        using aabb_t = axis_aligned_bounding_volume<cuboid3D<>>;

        // Only build box portals for cuboid volumes
        assert(volume.id() == volume_id::e_cuboid);
        // Need surfaces to wrap
        std::size_t n_surfaces{surfaces.size()};
        assert(n_surfaces != 0u);
        // Make sure data is consistent
        assert(n_surfaces == transforms.size() and
               n_surfaces == masks.total_size());

        // The surfaces container is prefilled with other surfaces
        dindex surfaces_offset = static_cast<dindex>(n_surfaces);

        // Fetch the position in the mask tuple for the rectangle portals
        constexpr auto rectangle_id{detector_t::masks::id::e_portal_rectangle2};

        // The material will be added in a later step
        constexpr auto no_material = surface_t::material_id::e_none;
        material_link_t material_link{
            no_material,
            detail::invalid_value<typename material_link_t::index_type>()};

        // Max distance in case of infinite bounds
        constexpr scalar max_shift{0.01f * std::numeric_limits<scalar>::max()};

        // The bounding boxes around the module surfaces
        std::vector<aabb_t> boxes;
        boxes.reserve(n_surfaces);

        for (const auto &sf : surfaces) {
            masks.template visit<bounding_box_creator>(
                sf.mask(), m_envelope, transforms[sf.transform()], boxes);
        }
        // Build an aabb in the global space around the surface aabbs
        aabb_t world_box{boxes, boxes.size(), m_envelope};

        // translation
        const point3_t center = world_box.template center<point3_t>();

        // The world box local frame is the global coordinate frame
        const point3_t box_min = world_box.template loc_min<point3_t>();
        const point3_t box_max = world_box.template loc_max<point3_t>();

        // Get the half lengths for the rectangle sides and translation
        const point3_t h_lengths = 0.5f * (box_max - box_min);
        const scalar h_x{math_ns::abs(h_lengths[0])};
        const scalar h_y{math_ns::abs(h_lengths[1])};
        const scalar h_z{math_ns::abs(h_lengths[2])};

        // Volume links for the portal descriptors and the masks
        const dindex volume_idx{volume.index()};
        const nav_link_t volume_link{detail::invalid_value<nav_link_t>()};

        // Construct portals in the...

        //
        // ... x-y plane
        //
        // Only one rectangle needed for both surfaces
        mask_link_t mask_link{rectangle_id,
                              masks.template size<rectangle_id>()};
        masks.template emplace_back<rectangle_id>(empty_context{}, volume_link,
                                                  h_x, h_y);

        // No rotation, but shift in z for both faces
        vector3_t shift{0.f, 0.f, std::isinf(h_z) ? max_shift : h_z};
        transforms.emplace_back(ctx, static_cast<vector3_t>(center + shift));
        transforms.emplace_back(ctx, static_cast<vector3_t>(center - shift));

        // Build the portal surfaces
        dindex trf_idx{transforms.size(ctx) - 2};
        surfaces.emplace_back(trf_idx, mask_link, material_link, volume_idx,
                              dindex_invalid, surface_id::e_portal);

        surfaces.emplace_back(++trf_idx, mask_link, material_link, volume_idx,
                              dindex_invalid, surface_id::e_portal);

        //
        // ... x-z plane
        //
        ++mask_link;
        masks.template emplace_back<rectangle_id>(empty_context{}, volume_link,
                                                  h_x, h_z);

        // Rotate by 90deg around x-axis, plus shift in y
        shift = {0.f, std::isinf(h_y) ? max_shift : h_y, 0.f};
        vector3_t new_x{1.f, 0.f, 0.f};
        vector3_t new_z{0.f, -1.f, 0.f};
        transforms.emplace_back(ctx, static_cast<vector3_t>(center + shift),
                                new_z, new_x);
        transforms.emplace_back(ctx, static_cast<vector3_t>(center - shift),
                                new_z, new_x);

        surfaces.emplace_back(++trf_idx, mask_link, material_link, volume_idx,
                              dindex_invalid, surface_id::e_portal);

        surfaces.emplace_back(++trf_idx, mask_link, material_link, volume_idx,
                              dindex_invalid, surface_id::e_portal);

        //
        // ... y-z plane
        //
        ++mask_link;
        masks.template emplace_back<rectangle_id>(empty_context{}, volume_link,
                                                  h_z, h_y);

        // Rotate by 90deg around y-axis, plus shift in x
        shift = {std::isinf(h_x) ? max_shift : h_x, 0.f, 0.f};
        new_x = {0.f, 0.f, -1.f};
        new_z = {1.f, 0.f, 0.f};
        transforms.emplace_back(ctx, static_cast<vector3_t>(center + shift),
                                new_z, new_x);
        transforms.emplace_back(ctx, static_cast<vector3_t>(center - shift),
                                new_z, new_x);

        surfaces.emplace_back(++trf_idx, mask_link, material_link, volume_idx,
                              dindex_invalid, surface_id::e_portal);

        surfaces.emplace_back(++trf_idx, mask_link, material_link, volume_idx,
                              dindex_invalid, surface_id::e_portal);

        return {surfaces_offset, static_cast<dindex>(surfaces.size())};
    }

    private:
    /// Portal envelope (min distance between portals and volume surfaces)
    scalar_t m_envelope;
};

}  // namespace detray
