/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/io/common/detail/type_traits.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/materials/material.hpp"
#include "detray/tools/detector_builder.hpp"
#include "detray/tools/material_builder.hpp"
#include "detray/tools/material_factory.hpp"

// System include(s)
#include <memory>
#include <string>
#include <vector>

namespace detray {

/// @brief Abstract base class for tracking geometry readers
template <class detector_t>
class homogeneous_material_reader : public reader_interface<detector_t> {

    using base_type = reader_interface<detector_t>;
    /// IO material ids do not need to coincide with the detector ids,
    /// they are shared with ACTS
    using material_type = io::detail::material_type;
    using scalar_type = typename detector_t::scalar_type;

    protected:
    /// Tag the reader as "homogeneous material"
    inline static const std::string tag = "homogeneous_material";

    public:
    /// Same constructors for this class as for base_type
    using base_type::base_type;

    protected:
    /// Deserialize the detector material @param det_mat_data from its IO
    /// payload
    static void deserialize(
        detector_builder<typename detector_t::metadata, volume_builder>&
            det_builder,
        typename detector_t::name_map& /*name_map*/,
        const detector_homogeneous_material_payload& det_mat_data) {
        using mat_types = typename detector_t::material_container::value_types;
        using material_id = typename detector_t::materials::id;

        // Deserialize the material volume by volume
        for (const auto& mv_data : det_mat_data.volumes) {
            // Decorate the current volume builder with material
            auto vm_builder = det_builder.template decorate<material_builder>(
                static_cast<dindex>(mv_data.index));

            // Add the material data to the factory
            auto mat_factory = std::make_shared<material_factory<detector_t>>();
            for (const auto& slab_data : mv_data.mat_slabs) {

                mat_factory->add_material(material_id::e_slab,
                                          deserialize(slab_data));
            }
            if constexpr (mat_types::n_types == 2u) {
                if (mv_data.mat_rods.has_value()) {
                    for (const auto& rod_data : *(mv_data.mat_rods)) {

                        mat_factory->add_material(material_id::e_rod,
                                                  deserialize(rod_data));
                    }
                }
            }

            // Add the material to the volume
            vm_builder->add_sensitives(mat_factory);
        }
    }

    /// @returns material data for a material factory from a slab io payload
    /// @param slab_data
    static material_data<scalar_type> deserialize(
        const material_slab_payload& slab_data) {

        return {static_cast<scalar_type>(slab_data.thickness),
                deserialize(slab_data.mat)};
    }

    /// @returns the material from its IO payload @param mat_data
    static auto deserialize(const material_payload& mat_data) {

        return material<scalar_type>{
            static_cast<scalar_type>(mat_data.params[0]),
            static_cast<scalar_type>(mat_data.params[1]),
            static_cast<scalar_type>(mat_data.params[2]),
            static_cast<scalar_type>(mat_data.params[3]),
            static_cast<scalar_type>(mat_data.params[4]),
            // The molar density is calculated on the fly
            static_cast<material_state>(mat_data.params[6])};
    }
};

}  // namespace detray
