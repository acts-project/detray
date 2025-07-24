/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/builders/homogeneous_material_factory.hpp"
#include "detray/builders/homogeneous_material_generator.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/builders/volume_builder_interface.hpp"
#include "detray/utils/log.hpp"

// System include(s)
#include <memory>
#include <stdexcept>
#include <vector>

namespace detray {

/// @brief Build a volume containing surfaces with material.
///
/// Decorator class to a volume builder that adds the material data to the
/// surfaces while building the volume.
template <typename detector_t>
class homogeneous_material_builder final : public volume_decorator<detector_t> {

    public:
    using material_id = typename detector_t::materials::id;
    using scalar_type = dscalar<typename detector_t::algebra_type>;

    /// @param vol_builder volume builder that should be decorated with material
    DETRAY_HOST
    explicit homogeneous_material_builder(
        std::unique_ptr<volume_builder_interface<detector_t>> vol_builder)
        : volume_decorator<detector_t>(std::move(vol_builder)) {}

    /// Overwrite, to add material in addition to surfaces (only if surfaces are
    /// present in the factory, otherwise only add material)
    /// @{
    DETRAY_HOST
    void add_surfaces(
        std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
        typename detector_t::geometry_context ctx = {}) override {

        DETRAY_DEBUG("homogeneous_material_builder::add_surface()");

        // If the factory also holds surface data, call base volume builder
        volume_decorator<detector_t>::add_surfaces(sf_factory, ctx);

        // Add material
        auto mat_factory =
            std::dynamic_pointer_cast<homogeneous_material_factory<detector_t>>(
                sf_factory);
        if (mat_factory) {
            DETRAY_DEBUG("-> found material factory: calling");
            (*mat_factory)(this->surfaces(), m_materials);
            return;
        }
        auto mat_generator = std::dynamic_pointer_cast<
            homogeneous_material_generator<detector_t>>(sf_factory);
        if (mat_generator) {
            DETRAY_DEBUG("-> found material generator: calling");
            (*mat_generator)(this->surfaces(), m_materials);
            return;
        }
    }
    /// @}

    /// Add the volume and the material to the detector @param det
    DETRAY_HOST
    auto build(detector_t &det, typename detector_t::geometry_context ctx = {})
        -> typename detector_t::volume_type * override {
        DETRAY_DEBUG("homogeneous_material_builder::build()");

        const auto &material = det.material_store();

        DETRAY_DEBUG("n_surfaces=" << this->surfaces().size());

        // Update the surface material links and shift them according to the
        // number of material slabs/rods that were in the detector previously
        for (auto &sf : this->surfaces()) {
            DETRAY_DEBUG("- sf=" << sf);
            if (sf.material().id() == material_id::e_slab) {
                dindex offset = material.template size<material_id::e_slab>();
                DETRAY_DEBUG("-> update material slab offset: " << offset);
                sf.update_material(offset);
            }
            if constexpr (detector_t::materials::template is_defined<
                              material_rod<scalar_type>>()) {
                if (sf.material().id() == material_id::e_rod) {
                    dindex offset =
                        material.template size<material_id::e_rod>();
                    DETRAY_DEBUG("-> update material rod offset: " << offset);
                    sf.update_material(
                        material.template size<material_id::e_rod>());
                }
            }
        }

        // Add material to the detector
        DETRAY_DEBUG("Appending "
                     << m_materials.template size<material_id::e_slab>()
                     << " slabs into detector materials");
        DETRAY_DEBUG("Appending "
                     << m_materials.template size<material_id::e_rod>()
                     << " rods into detector materials");
        det._materials.append(std::move(m_materials));
        m_materials.clear_all();

        // Call the underlying volume builder(s) and give the volume to the
        // next decorator
        return volume_decorator<detector_t>::build(det, ctx);
    }

    private:
    // Material container for this volume
    typename detector_t::material_container m_materials{};
};

}  // namespace detray
