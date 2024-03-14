/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/surface_factory_interface.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/geometry/shapes/cuboid3D.hpp"
#include "detray/utils/bounding_volume.hpp"

// System include(s)
#include <cassert>
#include <iostream>
#include <limits>

namespace detray {

/// @brief configuration for the cylinder portal generator
template <typename scalar_t>
struct cylinder_portal_config {
    /// Autofit the portals using the layer radius and the module surfaces
    bool m_autofit{true};
    /// Build inner cylinder portal (will use the same distance to the layer
    /// that was found for the outer cylinder portal)
    bool m_build_inner{true};
    /// Minimal envelope for the portals (used in autofitting)
    scalar_t m_envelope{100.f * unit<scalar_t>::um};
    /// FIxed outer radius during autofit
    scalar_t m_fixed_outer_r{0.f};
    /// FIxed length of the cylinder
    scalar_t m_fixed_z{0.f};
    /// Don't autofit the cylinder radii, build from explicit radii
    std::vector<dindex> m_volume_links{dindex_invalid, dindex_invalid,
                                       dindex_invalid, dindex_invalid};
    /// Don't autofit the cylinder radii, build from explicit radii
    std::vector<scalar_t> m_radii{};
    /// Don't autofit the cylinders length, build from explicit z extent
    std::vector<scalar_t> m_z_positions{};

    /// Setters
    /// @{
    constexpr cylinder_portal_config &autofit(const bool b) {
        m_autofit = b;
        return *this;
    }
    constexpr cylinder_portal_config &build_inner(const bool b) {
        m_build_inner = b;
        return *this;
    }
    constexpr cylinder_portal_config &envelope(const scalar_t e) {
        m_envelope = e;
        return *this;
    }
    constexpr cylinder_portal_config &fixed_outer_radius(const scalar_t r) {
        m_fixed_outer_r = r;
        return *this;
    }
    constexpr cylinder_portal_config &fixed_half_length(const scalar_t z) {
        m_fixed_z = z;
        return *this;
    }
    constexpr cylinder_portal_config &link_north(const dindex l) {
        m_volume_links[0] = l;
        return *this;
    }
    constexpr cylinder_portal_config &link_south(const dindex l) {
        m_volume_links[1] = l;
        return *this;
    }
    constexpr cylinder_portal_config &link_east(const dindex l) {
        m_volume_links[2] = l;
        return *this;
    }
    constexpr cylinder_portal_config &link_west(const dindex l) {
        m_volume_links[3] = l;
        return *this;
    }
    /// @}

    /// Getters
    /// @{
    constexpr scalar_t envelope() const { return m_envelope; }
    constexpr bool autofit() const { return m_autofit; }
    constexpr bool build_inner() const { return m_build_inner; }
    constexpr scalar_t fixed_outer_radius() const { return m_fixed_outer_r; }
    constexpr scalar_t fixed_half_length() const { return m_fixed_z; }
    constexpr const auto &volume_links() const { return m_volume_links; }
    /// @}
};

/// @brief Generates a portal box around a volume that already contains surfaces
///
/// @tparam detector_t the type of detector the volume belongs to.
template <typename detector_t>
class cylinder_portal_generator final
    : public surface_factory_interface<detector_t> {

    using scalar_t = typename detector_t::scalar_type;
    using point3_t = typename detector_t::point3;
    using vector3_t = typename detector_t::vector3;
    using transform3_t = typename detector_t::transform3;

    /// A functor to construct global bounding boxes around masks
    struct bounding_box_creator {

        using aabb_t = axis_aligned_bounding_volume<cuboid3D>;

        template <typename mask_group_t, typename index_t>
        DETRAY_HOST_DEVICE inline void operator()(
            const mask_group_t &mask_group, const index_t &index,
            const scalar_t envelope, const transform3_t &trf,
            std::vector<aabb_t> &boxes) const {
            // Local minimum bounding box
            aabb_t box{mask_group.at(index), boxes.size(), envelope};
            // Bounding box in global coordinates (might no longer be minimum)
            boxes.push_back(box.transform(trf));
        }
    };

    public:
    /// Save the boundaries of the cylinder after autofitting the portals
    struct boundaries {
        scalar_t inner_radius{0.f}, outer_radius{0.f}, lower_z{0.f},
            upper_z{0.f};
    };

    /// Use @param env as portal envelope
    DETRAY_HOST
    cylinder_portal_generator(const cylinder_portal_config<scalar_t> cfg)
        : m_cfg{cfg} {}

    /// @returns the number of rectangle portals this factory will produce
    DETRAY_HOST
    auto size() const -> dindex override {
        return m_cfg.autofit() ? 4u : static_cast<dindex>(m_portals.size());
    }

    DETRAY_HOST
    void clear() override {
        m_portals.clear();
        m_transforms.clear();
        m_cyl_bounds.clear();
        m_disc_bounds.clear();
    };

    DETRAY_HOST
    void push_back(surface_data<detector_t> &&) override { /*Do nothing*/
    }
    DETRAY_HOST
    auto push_back(std::vector<surface_data<detector_t>> &&)
        -> void override { /*Do nothing*/
    }

    /// @brief Add a single cylinder portal
    void add_cylinder_portal(dindex vol_link, const scalar_t r,
                             const scalar_t lower_z, const scalar_t upper_z) {
        // Call cylinder creation while ensuring that auto-construction is off
        add_cylinder_portal(vol_link, r, lower_z, upper_z, false);
    }

    /// @brief Add a single disc portal
    void add_disc_portal(dindex vol_link, const scalar_t inner_r,
                         const scalar_t outer_r, const scalar_t z) {
        // Call disc creation while ensuring that auto-construction is off
        add_disc_portal(vol_link, inner_r, outer_r, z, false);
    }

    /// @brief Access the volume boundaries after fitting
    const boundaries &volume_boundaries() { return m_boundaries; }

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
                    typename detector_t::surface_lookup_container &surfaces,
                    typename detector_t::transform_container &transforms,
                    typename detector_t::mask_container &masks,
                    typename detector_t::geometry_context ctx = {})
        -> dindex_range override {

        // Only build portals for cylinder volumes
        assert(volume.id() == volume_id::e_cylinder);
        const std::size_t n_surfaces{surfaces.size()};

        if (m_cfg.autofit()) {

            if (!m_transforms.empty() || !m_cyl_bounds.empty() ||
                !m_disc_bounds.empty() || !m_portals.empty()) {
                throw std::runtime_error(
                    "Factory already filled with data. Cannot do an automatic "
                    "portal fit!");
            }

            // Fill the data stores of the factory with portals automatically
            autofit_portals(ctx, surfaces, transforms, masks);
        }

        constexpr auto invalid_src_link{detail::invalid_value<std::uint64_t>()};

        // Build the portal surfaces
        for (auto &pt_desc : m_portals) {
            // Adjust the volume index
            pt_desc.set_volume(volume.index());
            // Update the mask link
            masks.template visit<detail::mask_index_update>(pt_desc.mask(),
                                                            pt_desc);
            // Update transform link
            pt_desc.update_transform(transforms.size(ctx));
            // Add to volume builder
            surfaces.push_back(pt_desc, invalid_src_link);
        }

        // Append the transforms
        transforms.insert(std::move(m_transforms), ctx);

        // Appendthe masks
        constexpr auto cyl_id{detector_t::masks::id::e_portal_cylinder2};
        constexpr auto disc_id{detector_t::masks::id::e_portal_ring2};

        for (const auto &cyl_bounds : m_cyl_bounds) {
            masks.template emplace_back<cyl_id>(
                empty_context{}, cyl_bounds.second, cyl_bounds.first);
        }

        for (const auto &disc_bounds : m_disc_bounds) {
            masks.template emplace_back<disc_id>(
                empty_context{}, disc_bounds.second, disc_bounds.first);
        }

        // Clear out all data
        clear();

        return {static_cast<dindex>(n_surfaces),
                static_cast<dindex>(surfaces.size())};
    }

    private:
    /// @brief Add a single cylinder portal to the internal data storage
    void add_cylinder_portal(dindex vol_link, const scalar_t r,
                             const scalar_t lower_z, const scalar_t upper_z,
                             bool do_autofit = true) {

        // Toggle autofit
        m_cfg.autofit(do_autofit);

        using surface_type = typename detector_t::surface_type;
        using mask_id = typename detector_t::masks::id;
        using mask_link_type = typename surface_type::mask_link;
        using material_id = typename detector_t::materials::id;
        using material_link_type = typename surface_type::material_link;

        const scalar_t min_z{math::min(lower_z, upper_z)};
        const scalar_t max_z{math::max(lower_z, upper_z)};

        // translation
        const point3_t tsl{0.f, 0.f, 0.f};

        // Add transform and mask data
        m_transforms.emplace_back(tsl);
        m_cyl_bounds.emplace_back(vol_link,
                                  std::vector<scalar_t>{r, min_z, max_z});

        // Add surface links
        mask_link_type mask_link{mask_id::e_portal_cylinder2,
                                 static_cast<dindex>(m_cyl_bounds.size() - 1u)};
        material_link_type material_link{material_id::e_none, dindex_invalid};

        m_portals.emplace_back(static_cast<dindex>(m_transforms.size() - 1u),
                               mask_link, material_link, dindex_invalid,
                               surface_id::e_portal);
    }

    /// @brief Add a single disc portal
    void add_disc_portal(dindex vol_link, const scalar_t inner_r,
                         const scalar_t outer_r, const scalar_t z,
                         bool do_autofit = true) {

        // Toggle autofit
        m_cfg.autofit(do_autofit);

        using surface_type = typename detector_t::surface_type;
        using mask_id = typename detector_t::masks::id;
        using mask_link_type = typename surface_type::mask_link;
        using material_id = typename detector_t::materials::id;
        using material_link_type = typename surface_type::material_link;

        const scalar_t min_r{math::min(inner_r, outer_r)};
        const scalar_t max_r{math::max(inner_r, outer_r)};

        // translation
        point3_t tsl{0.f, 0.f, z};

        // Add transform and mask data
        m_transforms.emplace_back(tsl);
        m_disc_bounds.emplace_back(vol_link,
                                   std::vector<scalar_t>{min_r, max_r});

        // Add surface links
        mask_link_type mask_link{mask_id::e_portal_ring2,
                                 static_cast<dindex>(m_cyl_bounds.size() - 1u)};
        material_link_type material_link{material_id::e_none, dindex_invalid};

        m_portals.emplace_back(static_cast<dindex>(m_transforms.size() - 1u),
                               mask_link, material_link, dindex_invalid,
                               surface_id::e_portal);
    }

    /// @brief Add a single disc portal
    template <typename surface_container_t, typename transform_container_t>
    scalar_t get_mean_radius(const surface_container_t &surfaces,
                             const transform_container_t &transforms) {

        double mean{0.};

        for (const auto &sf_desc : surfaces) {
            const auto &trf = transforms[sf_desc.transform()];
            mean += getter::perp(trf.translation());
        }

        return static_cast<scalar_t>(mean /
                                     static_cast<double>(surfaces.size()));
    }

    /// Automatically find the correct portal dimensions from fitting an
    /// axis aligned bounding box aroud the layer.
    void autofit_portals(
        typename detector_t::geometry_context ctx,
        typename detector_t::surface_lookup_container &surfaces,
        typename detector_t::transform_container &transforms,
        typename detector_t::mask_container &masks) {

        using aabb_t = axis_aligned_bounding_volume<cuboid3D>;

        // Need surfaces to wrap
        const std::size_t n_surfaces{surfaces.size()};
        assert(n_surfaces != 0u);
        // Make sure data is consistent
        assert(n_surfaces == transforms.size() and
               n_surfaces == masks.total_size());

        // The bounding boxes around the module surfaces
        std::vector<aabb_t> boxes;
        boxes.reserve(n_surfaces);

        for (const auto &sf : surfaces) {
            masks.template visit<bounding_box_creator>(
                sf.mask(), m_cfg.envelope(), transforms.at(sf.transform(), ctx),
                boxes);
        }

        // Build an aabb in the global space around the surface aabbs
        aabb_t world_box{boxes, boxes.size(), m_cfg.envelope()};

        // The world box local frame is the global coordinate frame
        const point3_t box_min = world_box.template loc_min<point3_t>();
        const point3_t box_max = world_box.template loc_max<point3_t>();

        // Get the half lengths for the cylinder height and disc translation
        const point3_t h_lengths = 0.5f * (box_max - box_min);
        const scalar h_x{math::abs(h_lengths[0])};
        const scalar h_y{math::abs(h_lengths[1])};
        const scalar h_z{math::abs(h_lengths[2])};

        const scalar_t outer_r_min{math::max(h_x, h_y)};
        const scalar_t mean_radius{get_mean_radius(surfaces, transforms)};

        scalar_t outer_r{outer_r_min};
        scalar_t inner_r{mean_radius - (outer_r - mean_radius)};
        scalar_t lower_z{-h_z};
        scalar_t upper_z{h_z};

        if (m_cfg.fixed_outer_radius() > 0.f) {
            outer_r = m_cfg.fixed_outer_radius();
        }
        if (m_cfg.fixed_half_length() > 0.f) {
            lower_z = -m_cfg.fixed_half_length();
            upper_z = m_cfg.fixed_half_length();
        }

        m_boundaries = {inner_r, outer_r, lower_z, upper_z};

        // If inner radius is 0, skip adding the inner cylinder
        if (m_cfg.build_inner()) {
            add_cylinder_portal(m_cfg.m_volume_links[1], inner_r, lower_z,
                                upper_z, true);
        }

        add_cylinder_portal(m_cfg.m_volume_links[0], outer_r, lower_z, upper_z,
                            true);

        add_disc_portal(m_cfg.m_volume_links[2], inner_r, outer_r, lower_z,
                        true);
        add_disc_portal(m_cfg.m_volume_links[3], inner_r, outer_r, upper_z,
                        true);
    }

    /// Portal generator configuration
    cylinder_portal_config<scalar_t> m_cfg;
    ///  Save the domensions of the volume after autofitting
    boundaries m_boundaries{};
    /// Portal descriptors
    std::vector<typename detector_t::surface_type> m_portals{};
    /// Transforms of surfaces
    dvector<typename detector_t::transform3> m_transforms{};
    /// Mask boundaries of surfaces
    std::vector<std::pair<dindex, std::vector<scalar_t>>> m_cyl_bounds{},
        m_disc_bounds{};
};

}  // namespace detray
