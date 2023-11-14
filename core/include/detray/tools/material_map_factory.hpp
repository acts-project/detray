/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/masks/unmasked.hpp"
#include "detray/materials/material.hpp"
#include "detray/tools/homogeneous_material_factory.hpp"
#include "detray/tools/surface_factory.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <vector>

namespace detray {

/// @brief Factory class for homogeneous material.
///
/// Uses a surface factory underneath the hood that handles the surface
/// construction. The material ID from the detector determines the type of
/// material that will be produced. The factory is filled from [a vector] of
/// @c material_data .
///
/// @tparam detector_t type of detector that contains the material
template <typename detector_t, typename index_t = dindex>
class material_map_factory final : public factory_decorator<detector_t> {

    using mask_id = typename detector_t::masks::id;
    using material_id = typename detector_t::materials::id;
    using index_type = index_t;

    using base_factory = factory_decorator<detector_t>;
    using placeholder_factory_t = surface_factory<detector_t, unmasked>;

    public:
    using scalar_type = typename detector_t::scalar_type;
    using data_type = material_data<scalar_type>;

    /// Factory with surfaces potentially already filled or empty placeholder
    /// that will not be used.
    DETRAY_HOST
    material_map_factory(
        std::unique_ptr<surface_factory_interface<detector_t>> sf_factory =
            std::make_unique<placeholder_factory_t>())
        : base_factory(std::move(sf_factory)) {}

    /// @returns the number of material instances that will be built by the
    /// factory
    DETRAY_HOST
    auto n_materials() const -> dindex {

        const dindex n_surfaces{static_cast<dindex>(m_links.size())};

        // Need exactly one material per surface
        assert(m_indices.empty() or (m_indices.size() == n_surfaces));
        assert(m_links.size() == n_surfaces);
        assert(m_materials.size() == n_surfaces);
        assert(m_thickness.size() == n_surfaces);

        return n_surfaces;
    }

    /// @returns the material links to the surfaces (counted for this volume)
    DETRAY_HOST
    auto links() const
        -> const std::map<std::size_t,
                          std::pair<material_id, std::vector<index_type>>> & {
        return m_links;
    }

    /// @returns the raw materials that are currently in the factory
    DETRAY_HOST
    auto materials() const
        -> const std::map<std::size_t, std::vector<material<scalar_type>>> & {
        return m_materials;
    }

    /// @returns the material thickness currently held by the factory
    DETRAY_HOST
    auto thickness() const
        -> const std::map<std::size_t, std::vector<scalar_type>> & {
        return m_thickness;
    }

    /// Add all necessary compontents to the factory for a material map
    ///
    /// @param mat_data contains vectors of material parameters
    /// @param indices local bin indices in the same order as the material data
    DETRAY_HOST
    void add_material(material_id id, data_type &&mat_data,
                      std::vector<std::size_t> n_bins,
                      std::vector<index_type> indices) {

        auto [sf_index, mat, thickness] = mat_data.get_data();

        m_links[sf_index] = std::make_pair(id, std::move(indices));
        m_indices.push_back(sf_index);
        m_n_bins[sf_index] = std::move(n_bins);
        m_materials[sf_index] = std::move(mat);
        m_thickness[sf_index] = std::move(thickness);
    }

    /// @brief Add data for material maps to the containers of a volume builder.
    ///
    /// This assumes that the surfaces have already been added (e.g. by this
    /// factories underlying surface factory).
    ///
    /// @note This does not override the pure virtual function from the surface
    /// factory interface, but presents an overload for the case when material
    /// should be added.
    ///
    /// @param surfaces surface container of the volume builder that should get
    ///                 decorated with material maps.
    /// @param materials map the local bin indices and their content to a
    ///                  volume local surface index.
    template <typename bin_data_t, std::size_t N>
    DETRAY_HOST auto operator()(
        typename detector_t::surface_lookup_container &surfaces,
        std::map<dindex, std::vector<bin_data_t>> &materials,
        std::map<dindex, std::array<std::size_t, N>> &n_bins) {

        using link_t = typename detector_t::surface_type::material_link;

        // Nothing left to do
        if (m_materials.empty() or m_thickness.empty()) {
            return;
        }

        // Add the material only to those surfaces that the data links against
        for (auto [sf_idx, sf] : detray::views::pick(surfaces, m_indices)) {

            assert(m_n_bins.at(sf_idx).size() == N);
            std::copy_n(m_n_bins.at(sf_idx).begin(), N, n_bins[sf_idx].begin());

            // Combine the material slab with its local bin index
            for (const auto [m_i, mat] :
                 detray::views::enumerate(m_materials.at(sf_idx))) {

                scalar_type t = m_thickness.at(sf_idx)[m_i];
                material_slab<scalar_type> mat_slab{mat, t};
                bin_data_t data{m_links[sf_idx].second.at(m_i), mat_slab};

                auto search = materials.find(sf_idx);
                if (search == materials.end()) {
                    materials[sf_idx] = std::vector<bin_data_t>{data};
                } else {
                    search->second.push_back(data);
                }
            }

            surfaces[sf_idx].material() =
                link_t{m_links[sf_idx].first, dindex_invalid};
        }
    }

    protected:
    /// Type and position(s) of the material in the detector material collection
    std::map<std::size_t, std::pair<material_id, std::vector<index_type>>>
        m_links{};
    /// Number of bins for the material grid axes
    std::map<std::size_t, std::vector<std::size_t>> m_n_bins{};
    /// Indices to link each material(vector) to its surface
    std::vector<std::size_t> m_indices{};
    /// Material thickness per surface
    std::map<std::size_t, std::vector<scalar_type>> m_thickness{};
    /// The pre-computed material to be wrapped in a slab or rod per surface
    std::map<std::size_t, std::vector<material<scalar_type>>> m_materials{};
};

}  // namespace detray
