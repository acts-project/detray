/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/homogeneous_material_builder.hpp"
#include "detray/builders/homogeneous_material_factory.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/io/common/detail/type_info.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/frontend/utils/type_traits.hpp"
#include "detray/materials/material.hpp"

// System include(s)
#include <memory>
#include <string>
#include <vector>

namespace detray::io {

template <typename detector_t, typename DIM>
class material_map_reader;

/// @brief Abstract base class for a homogeneous material reader.
template <class detector_t>
class homogeneous_material_reader : public reader_interface<detector_t> {

    using base_type = reader_interface<detector_t>;
    /// IO material ids do not need to coincide with the detector ids,
    /// they are shared with ACTS
    using material_type = io::material_id;
    using scalar_type = typename detector_t::scalar_type;

    friend class material_map_reader<detector_t,
                                     std::integral_constant<std::size_t, 2>>;

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

        using mat_id = typename detector_t::materials::id;

        // Deserialize the material volume by volume
        for (const auto& mv_data : det_mat_data.volumes) {
            // Decorate the current volume builder with material
            auto vm_builder = det_builder.template decorate<
                homogeneous_material_builder<detector_t>>(
                base_type::deserialize(mv_data.volume_link));

            // Add the material data to the factory
            auto mat_factory =
                std::make_shared<homogeneous_material_factory<detector_t>>();

            for (const auto& slab_data : mv_data.mat_slabs) {
                assert(slab_data.type == io::material_id::slab);

                const auto sf_link{
                    slab_data.index_in_coll.has_value()
                        ? slab_data.index_in_coll.value()
                        : detray::detail::invalid_value<std::size_t>()};

                mat_factory->add_material(mat_id::e_slab,
                                          deserialize(slab_data), sf_link);
            }
            if constexpr (detray::detail::has_material_rods_v<detector_t>) {
                if (mv_data.mat_rods.has_value()) {
                    for (const auto& rod_data : *(mv_data.mat_rods)) {
                        assert(rod_data.type == io::material_id::rod);

                        const auto sf_link{
                            rod_data.index_in_coll.has_value()
                                ? rod_data.index_in_coll.value()
                                : detray::detail::invalid_value<std::size_t>()};

                        mat_factory->add_material(
                            mat_id::e_rod, deserialize(rod_data), sf_link);
                    }
                }
            }

            // Add the material to the volume
            vm_builder->add_surfaces(mat_factory);
        }
    }

    /// @returns material data for a material factory from a slab io payload
    /// @param slab_data
    static material_data<scalar_type> deserialize(
        const material_slab_payload& slab_data) {

        return {static_cast<scalar_type>(slab_data.thickness),
                deserialize(slab_data.mat),
                base_type::deserialize(slab_data.surface)};
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

}  // namespace detray::io
