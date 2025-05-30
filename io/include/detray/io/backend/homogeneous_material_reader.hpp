/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/homogeneous_material_builder.hpp"
#include "detray/builders/homogeneous_material_factory.hpp"
#include "detray/core/concepts.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/io/backend/detail/basic_converter.hpp"
#include "detray/io/backend/detail/type_info.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/materials/material.hpp"

// System include(s)
#include <memory>
#include <string_view>
#include <vector>

namespace detray::io {

/// @brief Homogeneous material reader backend
///
/// Fills a @c detector_builder from a @c detector_homogeneous_material_payload
class homogeneous_material_reader {

    /// IO material ids do not need to coincide with the detector ids,
    /// they are shared with ACTS
    using material_type = io::material_id;

    public:
    /// Tag the reader as "homogeneous material"
    static constexpr std::string_view tag = "homogeneous_material";

    /// Payload type that the reader processes
    using payload_type = detector_homogeneous_material_payload;

    /// Convert the detector material @param det_mat_data from its IO
    /// payload
    template <class detector_t>
    static void from_payload(detector_builder<typename detector_t::metadata,
                                              volume_builder>& det_builder,
                             const payload_type& det_mat_data) {

        using scalar_t = dscalar<typename detector_t::algebra_type>;
        using mat_id = typename detector_t::materials::id;

        // Convert the material volume by volume
        for (const auto& mv_data : det_mat_data.volumes) {

            const auto vol_idx{
                detail::basic_converter::from_payload(mv_data.volume_link)};

            if (!det_builder.has_volume(vol_idx)) {
                std::stringstream err_stream;
                err_stream << "Volume " << vol_idx << ": "
                           << "Cannot build homogeneous material for volume "
                           << "(volume not registered in detector builder)";
                throw std::invalid_argument(err_stream.str());
            }

            // Decorate the current volume builder with material
            auto vm_builder = det_builder.template decorate<
                homogeneous_material_builder<detector_t>>(vol_idx);

            // Add the material data to the factory
            auto mat_factory =
                std::make_shared<homogeneous_material_factory<detector_t>>();

            for (const auto& slab_data : mv_data.mat_slabs) {
                assert(slab_data.type == io::material_id::slab);

                const auto sf_link{
                    slab_data.index_in_coll.has_value()
                        ? slab_data.index_in_coll.value()
                        : detray::detail::invalid_value<std::size_t>()};

                mat_factory->add_material(
                    mat_id::e_slab, from_payload<scalar_t>(slab_data), sf_link);
            }
            if constexpr (detray::concepts::has_material_rods<detector_t>) {
                if (mv_data.mat_rods.has_value()) {
                    for (const auto& rod_data : *(mv_data.mat_rods)) {
                        assert(rod_data.type == io::material_id::rod);

                        const auto sf_link{
                            rod_data.index_in_coll.has_value()
                                ? rod_data.index_in_coll.value()
                                : detray::detail::invalid_value<std::size_t>()};

                        mat_factory->add_material(
                            mat_id::e_rod, from_payload<scalar_t>(rod_data),
                            sf_link);
                    }
                }
            }

            // Add the material to the volume
            vm_builder->add_surfaces(mat_factory);
        }
    }

    /// @returns material data for a material factory from a slab io payload
    /// @param slab_data
    template <detray::concepts::scalar scalar_t>
    static material_data<scalar_t> from_payload(
        const material_slab_payload& slab_data) {

        return {static_cast<scalar_t>(slab_data.thickness),
                from_payload<scalar_t>(slab_data.mat),
                detail::basic_converter::from_payload(slab_data.surface)};
    }

    /// @returns the material from its IO payload @param mat_data
    template <detray::concepts::scalar scalar_t>
    static auto from_payload(const material_payload& mat_data) {

        return material<scalar_t>{
            static_cast<scalar_t>(mat_data.params[0]),
            static_cast<scalar_t>(mat_data.params[1]),
            static_cast<scalar_t>(mat_data.params[2]),
            static_cast<scalar_t>(mat_data.params[3]),
            static_cast<scalar_t>(mat_data.params[4]),
            // The molar density is calculated on the fly
            static_cast<material_state>(mat_data.params[6])};
    }
};

}  // namespace detray::io
