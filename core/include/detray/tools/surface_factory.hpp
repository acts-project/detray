/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/unmasked.hpp"
#include "detray/tools/surface_factory_interface.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <cassert>
#include <limits>
#include <memory>
#include <tuple>
#include <type_traits>
#include <vector>

namespace detray {

/// @brief Generates a number of surfaces for a volume (can be portals, passives
/// or sensitives) and fills them into the containers of a volume builder.
///
/// @tparam detector_t the type of detector the volume belongs to.
/// @tparam mask_shape_t the shape of the surface.
/// @tparam mask_id the concrete mask id that must be defined in the detector_t.
/// @tparam sf_id eihter portal, passive or sensitive.
template <typename detector_t, typename mask_shape_t,
          typename detector_t::mask_link::id_type mask_id, surface_id sf_id,
          typename volume_link_t = dindex>
class surface_factory final
    : public surface_factory_interface<detector_t>,
      public std::enable_shared_from_this<surface_factory<
          detector_t, mask_shape_t, mask_id, sf_id, volume_link_t>> {
    public:
    using scalar_t = typename detector_t::scalar_type;
    // Set individual volume link for portals, but only the mothervolume index
    // for other surfaces.
    using volume_link_collection =
        std::conditional_t<sf_id == surface_id::e_portal,
                           std::vector<volume_link_t>,
                           detray::views::single<volume_link_t>>;
    // ensures correct returning of shared pointers without pointless
    // construction and destructing of surface factories
    using smart_ptr_handler = std::enable_shared_from_this<surface_factory<
        detector_t, mask_shape_t, mask_id, sf_id, volume_link_t>>;

    /// shorthad for a colleciton of surface data that can be read by a surface
    /// factory
    // using surface_data_t = std::unique_ptr<
    //     surface_data_interface<detector_t, mask_shape_t, volume_link_t>>;
    using surface_data_t = surface_data<detector_t>;
    using sf_data_collection = std::vector<surface_data_t>;

    /// Empty factory.
    DETRAY_HOST
    surface_factory() = default;

    /// @returns the id for the surface type the factory handles (e.g. portals)
    DETRAY_HOST
    auto surface_type() const -> surface_id override { return sf_id; }

    /// @returns the current number of surfaces that will be built by this
    /// factory
    DETRAY_HOST
    auto size() const -> dindex override {
        check();
        return static_cast<dindex>(m_components.size());
    }

    /// @returns the mask boundaries currently held by the factory
    DETRAY_HOST
    auto components() const -> const std::vector<std::vector<scalar_t>> & {
        return m_components;
    }

    /// @returns the transforms currently held by the factory
    DETRAY_HOST
    auto transforms() const
        -> const std::vector<typename detector_t::transform3> & {
        return m_transforms;
    }

    /// @returns the volume link(s) currently held by the factory
    DETRAY_HOST
    const auto &volume_links() const { return m_volume_link; }

    /// Add all necessary compontents to the factory for a single surface
    DETRAY_HOST
    void push_back(surface_data_t &&sf_data) override {

        auto [trf, vlink, bounds] = sf_data.get_data();

        assert(bounds.size() == mask_shape_t::boundaries::e_size);

        if constexpr (sf_id != surface_id::e_portal) {
            if (size() == 0u) {
                *m_volume_link = vlink;
            } else {
                assert(*m_volume_link == vlink);
            }
        } else {
            m_volume_link.push_back(vlink);
        }

        m_transforms.push_back(trf);
        m_components.push_back(std::move(bounds));
    }

    /// Add all necessary compontents to the factory from bundled surface
    /// data in @param surface_data .
    DETRAY_HOST
    auto push_back(sf_data_collection &&surface_data)
        -> void /*std::shared_ptr<surface_factory<detector_t, mask_shape_t,
                   mask_id, sf_id, volume_link_t>>*/
        override {
        const auto n_surfaces{
            static_cast<dindex>(size() + surface_data.size())};

        m_transforms.reserve(n_surfaces);
        m_components.reserve(n_surfaces);

        if constexpr (sf_id == surface_id::e_portal) {
            m_volume_link.reserve(n_surfaces);
        }

        // Get per-surface data into detector level container layout
        for (auto &sf_data : surface_data) {
            push_back(std::move(sf_data));
        }

        // return smart_ptr_handler::shared_from_this();
    }

    /// Clear old data
    DETRAY_HOST
    auto clear() -> void {
        m_components.clear();
        m_transforms.clear();
        // cannot clear the single-view
        if constexpr (sf_id == surface_id::e_portal) {
            m_volume_link.clear();
        } else {
            *m_volume_link = volume_link_t{};
        }
    }

    /// Generate the surfaces and add them to given data collections.
    ///
    /// @param volume the volume they will be added to in the detector.
    /// @param surfaces the resulting surface objects.
    /// @param transforms the transforms of the surfaces.
    /// @param masks the masks of the surfaces (all of the same shape).
    /// @param ctx the geometry context.
    DETRAY_HOST
    auto operator()(const typename detector_t::volume_type &volume,
                    typename detector_t::surface_container_t &surfaces,
                    typename detector_t::transform_container &transforms,
                    typename detector_t::mask_container &masks,
                    typename detector_t::geometry_context ctx = {}) const
        -> dindex_range override {
        using surface_t = typename detector_t::surface_type;
        using mask_link_t = typename surface_t::mask_link;
        using material_link_t = typename surface_t::material_link;

        // The material will be added in a later step
        constexpr auto no_material = surface_t::material_id::e_none;
        // In case the surfaces container is prefilled with other surfaces
        dindex surfaces_offset = static_cast<dindex>(surfaces.size());

        // Nothing to construct
        if (size() == 0u) {
            return {surfaces_offset, surfaces_offset};
        }

        for (const auto [idx, comp] : detray::views::enumerate(m_components)) {

            // Add transform
            transforms.push_back(m_transforms[idx], ctx);

            if constexpr (std::is_same_v<mask_shape_t, unmasked>) {
                masks.template emplace_back<mask_id>(
                    empty_context{}, m_volume_link[idx],
                    std::numeric_limits<scalar_t>::infinity());
            } else {
                masks.template emplace_back<mask_id>(empty_context{}, comp,
                                                     m_volume_link[idx]);
            }

            // Add surface with all links set (relative to the given containers)
            mask_link_t mask_link{mask_id, masks.template size<mask_id>() - 1u};
            material_link_t material_link{no_material, dindex_invalid};
            surfaces.emplace_back(transforms.size(ctx) - 1u, mask_link,
                                  material_link, volume.index(), dindex_invalid,
                                  sf_id);
        }

        return {surfaces_offset, static_cast<dindex>(surfaces.size())};
    }

    private:
    /// Check that the containers have the same size
    DETRAY_HOST
    void check() const {
        // This should not happend (need same number and ordering of data)
        assert(m_components.size() == m_transforms.size());
        if constexpr (sf_id == surface_id::e_portal) {
            assert(m_components.size() == m_volume_link.size());
        } else {
            assert(m_volume_link.size() == 1);
        }
    }

    std::vector<std::vector<scalar_t>> m_components{};
    std::vector<typename detector_t::transform3> m_transforms{};
    volume_link_collection m_volume_link{};
};

}  // namespace detray
