/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/materials/material.hpp"
#include "detray/tools/surface_factory.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <cassert>
#include <memory>
#include <tuple>
#include <type_traits>
#include <vector>

namespace detray {

/// @brief Bind components for material together.
template <typename scalar_t>
class material_data {
    public:
    /// Construct from a predefined material
    ///
    /// @param mat predefined material, see
    ///            'detray/materials/predefined_materials.hpp'
    /// @param thickness of the material slab/rod
    DETRAY_HOST
    constexpr material_data(const scalar_t thickness,
                            const material<scalar_t> &mat)
        : m_mat{mat}, m_thickness{thickness} {}

    /// Construct from all parameters:
    ///
    /// @param thickness of the material slab/rod
    /// @param material_paramters0 x0 is the radiation length
    /// @param material_paramters1 l0 is the nuclear interaction length
    /// @param material_paramters2 ar is the relative atomic mass
    /// @param material_paramters3 z is the nuclear charge number
    /// @param material_paramters4 molarRho is the molar density
    /// @param state of the material (liquid, solid etc.)
    DETRAY_HOST
    constexpr material_data(
        const scalar_t thickness,
        const std::vector<scalar_t> &material_paramters,
        const material_state state = material_state::e_solid)
        : m_mat{material_paramters[0], material_paramters[1],
                material_paramters[2], material_paramters[3],
                material_paramters[4], state},
          m_thickness{thickness} {}

    /// @returns tuple based access to the contained material data.
    DETRAY_HOST
    constexpr auto get_data() -> std::tuple<material<scalar_t> &, scalar_t &> {
        return std::tie(m_mat, m_thickness);
    }

    private:
    /// The material parametrization
    material<scalar_t> m_mat{};
    /// Thickness/radius of the material slab/rod
    scalar_t m_thickness{0.f};
};

/// @brief Factory class for homogeneous material.
///
/// Uses a surface factory underneath the hood that handles the surface
/// construction. The material ID from the detector determines the type of
/// material that will be produced. The factory is filled from [a vector] of
/// @c material_data .
///
/// @tparam detector_t type of detector that contains the material
template <typename detector_t>
class material_factory final : public factory_decorator<detector_t> {

    using base_factory = factory_decorator<detector_t>;

    public:
    using material_id = typename detector_t::materials::id;
    using scalar_type = typename detector_t::scalar_type;

    /// Factory with surfaces potentially already filled.
    DETRAY_HOST
    material_factory(std::unique_ptr<surface_factory_interface<detector_t>>
                         sf_factory = nullptr)
        : base_factory(std::move(sf_factory)) {}

    /// @returns the number of material that will be built by the factory
    DETRAY_HOST
    auto size() const -> dindex {

        const dindex n_surfaces{static_cast<dindex>(m_links.size())};

        // Need exactly one material per surface
        assert(m_materials.size() == n_surfaces);
        assert(m_thickness.size() == n_surfaces);

        return n_surfaces;
    }

    /// @returns the material links to the surface (counted for this volume)
    DETRAY_HOST
    auto links() const -> const std::vector<material_id> & { return m_links; }

    /// @returns the raw materials that are currently in the factory
    DETRAY_HOST
    auto materials() const -> const std::vector<material<scalar_type>> & {
        return m_materials;
    }

    /// @returns the material thickness currently held by the factory
    DETRAY_HOST
    auto thickness() const -> const std::vector<scalar_type> & {
        return m_thickness;
    }

    /// Add all necessary compontents to the factory for a single material slab
    /// or rod (determined by the @param id)
    DETRAY_HOST
    void add_material(material_id id, material_data<scalar_type> &&mat_data) {

        auto [mat, thickness] = mat_data.get_data();

        m_links.push_back(id);
        m_materials.push_back(mat);
        m_thickness.push_back(thickness);
    }

    /// Add all necessary compontents to the factory for multiple material slabs
    /// or rods (determined by the @param id)
    DETRAY_HOST
    void add_material(material_id id,
                      std::vector<material_data<scalar_type>> &&mat_data_vec) {
        // Read the material containers
        m_links.reserve(m_links.size() + mat_data_vec.size());
        m_materials.reserve(m_materials.size() + mat_data_vec.size());
        m_thickness.reserve(m_thickness.size() + mat_data_vec.size());

        // Add the material components
        for (auto &mat_data : mat_data_vec) {
            this->add_material(id, std::move(mat_data));
        }
    }

    /// @brief Add material to the containers of a volume builder.
    ///
    /// This assumes that the surfaces have already been added (e.g. by this
    /// factories underlying surface factory). The material from this factory is
    /// added to the corresponding number of surfaces at the end of the calling
    /// volume builder.
    ///
    /// @note This does not override the pure virtual function from the surface
    /// factory interface, but presents an overload for the case when material
    /// should be added.
    ///
    /// @param surfaces surface store of the volume builder that should get
    ///                 decorated with material.
    /// @param material material store of the volume builder that the new
    ///                 materials get added to.
    DETRAY_HOST
    auto operator()(typename detector_t::surface_container_t &surfaces,
                    typename detector_t::material_container &materials) const {

        using link_t = typename detector_t::surface_type::material_link;

        // Check that the surfaces were set up correctly
        const dindex n_material{this->size()};
        assert(surfaces.size() >= n_material);

        // Range of surfaces for which to add material
        dindex_range r{static_cast<dindex>(surfaces.size() - n_material),
                       static_cast<dindex>(surfaces.size())};
        std::size_t sf_idx{0u};
        for (auto &sf : detray::ranges::subrange(surfaces, r)) {

            const material<scalar_type> &mat = m_materials.at(sf_idx);
            scalar_type t = m_thickness.at(sf_idx);

            dindex n{0u};
            if (m_links.at(sf_idx) == material_id::e_slab) {
                n = materials.template size<material_id::e_slab>();
                materials.template emplace_back<material_id::e_slab>({}, mat,
                                                                     t);
            } else {
                n = materials.template size<material_id::e_rod>();
                materials.template emplace_back<material_id::e_rod>({}, mat, t);
            }

            // Set the initial surface material link (will be updated when
            // added to the detector)
            sf.material() = link_t{m_links[sf_idx], n};
            ++sf_idx;
        }
    }

    protected:
    /// Material links of surfaces (currently only type ids)
    std::vector<material_id> m_links{};
    /// Material thickness
    std::vector<scalar_type> m_thickness{};
    /// The pre-computed material to be wrapped in a slab or rod
    std::vector<material<scalar_type>> m_materials{};
};

}  // namespace detray
