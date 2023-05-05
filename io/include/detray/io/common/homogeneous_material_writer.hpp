/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/detail/utils.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/io/common/writer_interface.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <string>
#include <string_view>

namespace detray {

/// @brief Abstract base class for simple material description writers
template <class detector_t>
class homogeneous_material_writer : public writer_interface<detector_t> {

    using base_type = writer_interface<detector_t>;

    protected:
    /// Tag the writer as "homogeneous_material"
    inline static const std::string tag = "homogeneous_material";

    public:
    /// Same constructors for this class as for base_type
    using base_type::base_type;

    protected:
    /// Serialize the header information into its payload
    static homogeneous_material_header_payload write_header(
        const detector_t& det, const std::string_view det_name) {
        using mat_types = typename detector_t::material_container::value_types;

        homogeneous_material_header_payload header_data;

        header_data.version = detail::get_detray_version();
        header_data.detector = det_name;
        header_data.tag = tag;
        header_data.date = detail::get_current_date();

        const auto& materials = det.material_store();
        header_data.n_slabs = materials.template size<mat_types::to_id(0u)>();
        header_data.n_rods = 0u;
        if constexpr (mat_types::n_types == 2u) {
            // The compiler looks at this code, even if the number of material
            // types is one. Therefore, "mat_types::n_types - 1" is safer to use
            header_data.n_rods =
                materials
                    .template size<mat_types::to_id(mat_types::n_types - 1)>();
        }

        return header_data;
    }

    /// Serialize the material description of a detector @param det into its io
    /// payload
    static detector_homogeneous_material_payload serialize(
        const detector_t& det) {
        using mat_types = typename detector_t::material_container::value_types;

        detector_homogeneous_material_payload dm_data;

        const auto& materials = det.material_store();

        // Serialize material slabs and rods
        for (const auto [idx, mat] : detray::views::enumerate(
                 materials.template get<mat_types::to_id(0u)>())) {
            dm_data.mat_slabs.push_back(serialize(mat, idx));
        }
        if constexpr (mat_types::n_types == 2u) {
            dm_data.mat_rods = {};
            for (const auto [idx, mat] : detray::views::enumerate(
                     materials.template get<mat_types::to_id(
                         mat_types::n_types - 1)>())) {
                dm_data.mat_rods->push_back(serialize(mat, idx));
            }
        }

        return dm_data;
    }

    /// Serialize surface material @param mat into its io payload
    static material_payload serialize(
        const material<typename detector_t::scalar_type>& mat) {
        material_payload mat_data;

        mat_data.params = {mat.X0(),
                           mat.L0(),
                           mat.Ar(),
                           mat.Z(),
                           mat.mass_density(),
                           mat.molar_density(),
                           static_cast<real_io>(mat.state())};
        return mat_data;
    }

    /// Serialize a surface material slab @param mat_slab into its io payload
    static material_slab_payload serialize(
        const material_slab<typename detector_t::scalar_type>& mat_slab,
        std::size_t idx) {
        material_slab_payload mat_data;

        mat_data.index = idx;
        mat_data.thickness = mat_slab.thickness();
        mat_data.mat = serialize(mat_slab.get_material());

        return mat_data;
    }

    /// Serialize a line material rod @param mat_rod into its io payload
    static material_slab_payload serialize(
        const material_rod<typename detector_t::scalar_type>& mat_rod,
        std::size_t idx) {
        material_slab_payload mat_data;

        mat_data.index = idx;
        mat_data.thickness = mat_rod.radius();
        mat_data.mat = serialize(mat_rod.get_material());

        return mat_data;
    }
};

}  // namespace detray
