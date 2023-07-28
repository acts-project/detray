/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/common/payloads.hpp"
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

        header_data.common = base_type::serialize(det_name, tag);

        const auto& materials = det.material_store();

        header_data.sub_header.emplace();
        auto& mat_sub_header = header_data.sub_header.value();
        mat_sub_header.n_slabs =
            materials.template size<mat_types::to_id(0u)>();
        mat_sub_header.n_rods = 0u;
        if constexpr (mat_types::n_types == 2u) {
            // The compiler looks at this code, even if the number of material
            // types is one. Therefore, "mat_types::n_types - 1" is safer to use
            mat_sub_header.n_rods =
                materials
                    .template size<mat_types::to_id(mat_types::n_types - 1)>();
        }

        return header_data;
    }

    /// Serialize the material description of a detector @param det into its io
    /// payload
    static detector_homogeneous_material_payload serialize(
        const detector_t& det, const typename detector_t::name_map&) {
        detector_homogeneous_material_payload dm_data;
        dm_data.volumes.reserve((det.volumes().size()));

        for (const auto& vol : det.volumes()) {
            dm_data.volumes.push_back(serialize(vol, det));
        }

        return dm_data;
    }

    /// Serialize the material description of a volume @param vol_desc into its
    /// io payload
    static material_volume_payload serialize(
        const typename detector_t::volume_type& vol_desc,
        const detector_t& det) {
        using material_type = material_slab_payload::material_type;

        material_volume_payload mv_data;
        mv_data.index = vol_desc.index();

        // Find all surfaces that belong to the volume
        for (const auto& sf_desc : det.surface_lookup()) {
            if (sf_desc.volume() != vol_desc.index()) {
                continue;
            }
            // Serialize material slabs and rods
            const auto sf = surface{det, sf_desc};
            const material_slab_payload mslp =
                sf.template visit_material<get_material_payload>();

            if (mslp.type == material_type::slab) {
                mv_data.mat_slabs.push_back(mslp);
            } else if (mslp.type == material_type::rod) {
                if (not mv_data.mat_rods.has_value()) {
                    mv_data.mat_rods.emplace();
                }
                mv_data.mat_rods->push_back(mslp);
            } else {
                throw std::runtime_error(
                    "Material could not be matched to payload (found type " +
                    std::to_string(static_cast<int>(mslp.type)) + ")");
            }
        }

        return mv_data;
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

        mat_data.type = material_slab_payload::material_type::slab;
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

        mat_data.type = material_slab_payload::material_type::rod;
        mat_data.index = idx;
        mat_data.thickness = mat_rod.radius();
        mat_data.mat = serialize(mat_rod.get_material());

        return mat_data;
    }

    private:
    /// Retrieve @c material_slab_payload from a material store element
    struct get_material_payload {
        template <typename material_group_t, typename index_t>
        inline auto operator()(const material_group_t& material_group,
                               const index_t& index) const {
            return homogeneous_material_writer<detector_t>::serialize(
                material_group[index], index);
        }
    };
};

}  // namespace detray
