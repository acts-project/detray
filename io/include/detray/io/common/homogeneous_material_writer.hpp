/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/io/common/detail/type_info.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/frontend/utils/type_traits.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/utils/type_list.hpp"

// System include(s)
#include <string>
#include <string_view>

namespace detray::io {

template <typename detector_t>
class material_map_writer;

/// @brief Abstract base class for simple material description writers
template <class detector_t>
class homogeneous_material_writer : public writer_interface<detector_t> {

    using base_type = writer_interface<detector_t>;
    using scalar_t = typename detector_t::scalar_type;

    friend class material_map_writer<detector_t>;

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

        homogeneous_material_header_payload header_data;

        header_data.common = base_type::serialize(det_name, tag);

        const auto& materials = det.material_store();

        header_data.sub_header.emplace();
        auto& mat_sub_header = header_data.sub_header.value();
        if constexpr (detray::detail::has_material_slabs_v<detector_t>) {
            mat_sub_header.n_slabs =
                materials.template size<detector_t::materials::id::e_slab>();
        }
        mat_sub_header.n_rods = 0u;
        if constexpr (detray::detail::has_material_rods_v<detector_t>) {
            mat_sub_header.n_rods =
                materials.template size<detector_t::materials::id::e_rod>();
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
        using material_type = material_slab_payload::mat_type;

        material_volume_payload mv_data;
        mv_data.volume_link = base_type::serialize(vol_desc.index());

        // Return early if the stores for homogeneous materials are empty
        using mat_id = typename detector_t::materials::id;

        // If this reader is called, the detector has at least material slabs
        if (det.material_store().template empty<mat_id::e_slab>()) {
            // Check for material rods that are present in e.g. wire chambers
            if constexpr (detray::detail::has_material_rods_v<detector_t>) {
                if (det.material_store().template empty<mat_id::e_rod>()) {
                    return mv_data;
                }
            } else {
                return mv_data;
            }
        }

        // Find all surfaces that belong to the volume and count them
        std::size_t sf_idx{0u}, slab_idx{0u}, rod_idx{0u};
        auto vol = detector_volume{det, vol_desc};
        for (const auto& sf_desc : vol.surfaces()) {

            // Serialize material slabs and rods
            const auto sf = surface{det, sf_desc};
            material_slab_payload mslp =
                sf.template visit_material<get_material_payload>(sf_idx);

            if (mslp.type == material_type::slab) {
                mslp.index_in_coll = slab_idx++;
                mv_data.mat_slabs.push_back(mslp);
            } else if (mslp.type == material_type::rod) {
                if (not mv_data.mat_rods.has_value()) {
                    mv_data.mat_rods.emplace();
                }
                mslp.index_in_coll = rod_idx++;
                mv_data.mat_rods->push_back(mslp);
            } else { /* material maps are handled by another writer */
            }
            ++sf_idx;
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
    template <template <typename> class material_t>
    static material_slab_payload serialize(const material_t<scalar_t>& mat,
                                           std::size_t sf_idx) {
        material_slab_payload mat_data;

        mat_data.type = io::detail::get_id<material_t<scalar_t>>();
        mat_data.surface = base_type::serialize(sf_idx);
        mat_data.thickness = mat.thickness();
        mat_data.mat = serialize(mat.get_material());

        return mat_data;
    }

    private:
    /// Retrieve @c material_slab_payload from a material store element
    struct get_material_payload {
        template <typename material_group_t, typename index_t>
        inline auto operator()(const material_group_t& material_group,
                               const index_t& index,
                               [[maybe_unused]] std::size_t sf_index) const {
            using material_t = typename material_group_t::value_type;

            constexpr bool is_slab =
                std::is_same_v<material_t, material_slab<scalar_t>>;
            constexpr bool is_rod =
                std::is_same_v<material_t, material_rod<scalar_t>>;

            if constexpr (is_slab or is_rod) {
                return homogeneous_material_writer<detector_t>::serialize(
                    material_group[index], sf_index);
            } else {
                return material_slab_payload{};
            }
        }
    };
};

}  // namespace detray::io
