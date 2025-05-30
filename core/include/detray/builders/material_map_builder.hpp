/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/builders/bin_fillers.hpp"
#include "detray/builders/material_map_factory.hpp"
#include "detray/builders/material_map_generator.hpp"
#include "detray/builders/surface_factory_interface.hpp"
#include "detray/builders/volume_builder_interface.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/materials/material_map.hpp"

// System include(s)
#include <cassert>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

namespace detray {

namespace detail {

struct material_coll_size;

template <typename materials_t>
struct add_sf_material_map;

}  // namespace detail

/// @brief Build the material maps for a given volume.
///
/// Decorator class to a volume builder that adds material maps to either
/// surfaces or volumes
template <typename detector_t, std::size_t DIM = 2u,
          typename mat_map_factory_t =
              material_grid_factory<typename detector_t::algebra_type>>
class material_map_builder final : public volume_decorator<detector_t> {
    using materials_t = typename detector_t::materials;

    public:
    using scalar_type = dscalar<typename detector_t::algebra_type>;
    using detector_type = detector_t;
    using value_type = material_slab<scalar_type>;
    using bin_data_type =
        typename fill_by_bin::template bin_data<DIM, value_type>;

    /// @param vol_builder volume builder that should be decorated with material
    /// maps
    DETRAY_HOST
    explicit material_map_builder(
        std::unique_ptr<volume_builder_interface<detector_t>> vol_builder)
        : volume_decorator<detector_t>(std::move(vol_builder)) {}

    /// @returns the raw materials and their local bin indices that are
    /// currently staged in the builder
    DETRAY_HOST
    auto data() const -> const std::map<dindex, std::vector<bin_data_type>>& {
        return m_bin_data;
    }

    /// Overwrite, to add material maps in addition to surfaces
    /// @{
    DETRAY_HOST
    void add_surfaces(
        std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
        typename detector_t::geometry_context ctx = {}) override {

        // If the factory also holds surface data, call base volume builder
        volume_decorator<detector_t>::add_surfaces(sf_factory, ctx);

        // Add material and bin data
        auto mat_factory = std::dynamic_pointer_cast<
            material_map_factory<detector_t, axis::multi_bin<DIM>>>(sf_factory);
        if (mat_factory) {
            (*mat_factory)(this->surfaces(), m_bin_data, m_n_bins,
                           m_axis_spans);
        }
        auto mat_generator =
            std::dynamic_pointer_cast<material_map_generator<detector_t>>(
                sf_factory);
        if (mat_generator) {
            (*mat_generator)(this->surfaces(), this->masks(), m_bin_data,
                             m_n_bins);
            return;
        }
    }
    /// @}

    /// Add the volume and the material maps to the detector @param det
    DETRAY_HOST
    auto build(detector_t& det, typename detector_t::geometry_context ctx = {})
        -> typename detector_t::volume_type* override {

        // Ensure the material links are correct BEFORE the surfaces are built
        // and potentially added to an acceleration data structure
        update_material_links(det);

        // Construct the surfaces
        typename detector_t::volume_type* vol_ptr =
            volume_decorator<detector_t>::build(det, ctx);

        // Build the material grid for every surface that has material
        dindex sf_idx = 0u;
        auto vol = tracking_volume{det, this->vol_index()};
        for (const auto& sf_desc : vol.surfaces()) {

            if (!surface_has_map(sf_idx)) {
                sf_idx++;
                continue;
            }

            // Construct and append the material map for a given surface shape
            darray<std::vector<scalar_type>, DIM> axis_spans{};
            auto axis_spans_itr = m_axis_spans.find(sf_idx);
            if (axis_spans_itr != m_axis_spans.end()) {
                axis_spans = axis_spans_itr->second;
            }

            // Construct and append the material map for a given surface shape
            auto sf = geometry::surface{det, sf_desc};
            [[maybe_unused]] auto [mat_id, mat_idx] = sf.template visit_mask<
                detail::add_sf_material_map<materials_t>>(
                m_factory, m_bin_data.at(sf_idx), m_n_bins.at(sf_idx),
                axis_spans, det._materials);

            // Make sure the linking was precomputed correctly
            assert(mat_id == sf_desc.material().id());
            assert(mat_idx == sf_desc.material().index());
            sf_idx++;
        }

        // Give the volume to the next decorator
        return vol_ptr;
    }

    private:
    /// Check whether a surface with a given index @param sf_idx should receive
    /// material from this builder
    bool surface_has_map(const dindex sf_idx) const {
        return m_bin_data.contains(sf_idx);
    }

    /// Set the correct global surface material link for this detector
    void update_material_links(const detector_t& det) {

        // The total number of surfaces that will be built by this builder
        const dindex n_surfaces{static_cast<dindex>(this->surfaces().size())};

        // Count the number of maps per material type to get the correct indices
        std::map<typename materials_t::id, dindex> mat_type_count;

        for (dindex sf_idx = 0u; sf_idx < n_surfaces; ++sf_idx) {
            if (!surface_has_map(sf_idx)) {
                continue;
            }

            auto& sf_desc = this->surfaces().at(sf_idx);

            const auto id{sf_desc.material().id()};
            if (!mat_type_count.contains(id)) {
                mat_type_count.emplace(id, 0u);
            } else {
                mat_type_count.at(id)++;
            }

            sf_desc.material().set_index(mat_type_count.at(id));
        }

        // Current sizes of the material stores to get global index offsets
        std::map<std::size_t, dindex> size_map;
        det._materials.template apply<detail::material_coll_size>(
            size_map, std::make_index_sequence<materials_t::n_types>{});

        // Update the counts with the detector offset
        for (dindex sf_idx = 0u; sf_idx < n_surfaces; ++sf_idx) {
            if (!surface_has_map(sf_idx)) {
                continue;
            }
            auto& sf_desc = this->surfaces().at(sf_idx);

            auto coll_idx{static_cast<std::size_t>(sf_desc.material().id())};
            sf_desc.material() += size_map.at(coll_idx);
        }
    }

    /// The surface this material map belongs to (index is volume local)
    std::map<dindex, std::vector<bin_data_type>> m_bin_data;
    /// Number of bins for the material grid axes
    std::map<dindex, darray<std::size_t, DIM>> m_n_bins{};
    /// The Axis spans for the material grid axes
    std::map<dindex, darray<std::vector<scalar_type>, DIM>> m_axis_spans{};
    /// Helper to generate empty grids
    mat_map_factory_t m_factory{};
};

namespace detail {

/// A functor to obtain the material collection sizes for every collection
struct material_coll_size {

    template <typename... coll_ts, std::size_t... I>
    DETRAY_HOST inline void operator()(std::map<std::size_t, dindex>& size_map,
                                       std::index_sequence<I...> /*seq*/,
                                       const coll_ts&... coll) const {
        (size_map.emplace(I, static_cast<dindex>(coll.size())), ...);
    }
};

/// A functor to add a material map to a surface
template <typename materials_t>
struct add_sf_material_map {

    template <typename mask_coll_t, typename index_range_t,
              typename mat_factory_t, typename bin_data_t, std::size_t DIM,
              typename material_store_t, concepts::scalar scalar_t>
    DETRAY_HOST inline std::pair<typename materials_t::id, dindex> operator()(
        [[maybe_unused]] const mask_coll_t& mask_coll,
        [[maybe_unused]] const index_range_t& index,
        [[maybe_unused]] const mat_factory_t& mat_factory,
        [[maybe_unused]] std::vector<bin_data_t>& bin_data,
        [[maybe_unused]] const darray<std::size_t, DIM>& n_bins,
        [[maybe_unused]] const darray<std::vector<scalar_t>, DIM>& axis_spans,
        [[maybe_unused]] material_store_t& mat_store) const {

        using mask_t = typename mask_coll_t::value_type;
        using mask_shape_t = typename mask_t::shape;

        constexpr bool is_line{
            std::is_same_v<mask_shape_t, detray::line_square> ||
            std::is_same_v<mask_shape_t, detray::line_circular>};

        // No material maps for line surfaces
        if constexpr (!is_line && mask_shape_t::dim == DIM) {
            // Map a grid onto the surface mask (the boundaries are taken from
            // the @c axis_spans variable, if it is not empty)
            mask_t sf_mask = {};
            if constexpr (concepts::interval<index_range_t>) {
                using index_t = typename index_range_t::index_type;

                // Find the true surface extent over all masks
                sf_mask = mask_coll.at(index.lower());

                if (index.size() > 1u) {
                    const index_range_t other_masks{
                        index.lower() + 1u,
                        static_cast<index_t>(index.size() - 1u)};

                    // Merge sub-masks
                    for (const auto& sub_mask :
                         detray::ranges::subrange(mask_coll, other_masks)) {
                        sf_mask = sf_mask + sub_mask;
                    }
                }
            } else {
                sf_mask = mask_coll.at(index);
            }

            auto mat_grid =
                mat_factory.new_grid(sf_mask, n_bins, {}, {}, axis_spans);

            // The detector only knows the non-owning grid types
            using non_owning_t =
                typename decltype(mat_grid)::template type<false>;
            static_assert(materials_t::template is_defined<non_owning_t>());

            // Add the material slabs to the grid
            for (const auto& bin : bin_data) {
                mat_grid.template populate<replace<>>(bin.local_bin_idx,
                                                      bin.single_element);
            }

            // Add the material grid to the detector
            constexpr auto gid{materials_t::template get_id<non_owning_t>()};
            mat_store.template push_back<gid>(mat_grid);

            // Return the index of the new material map
            return {gid,
                    static_cast<dindex>(mat_store.template size<gid>() - 1u)};
        } else {
            return {materials_t::id::e_none, dindex_invalid};
        }
    }
};

}  // namespace detail

}  // namespace detray
