/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/detail/grid_writer.hpp"
#include "detray/io/common/detail/type_traits.hpp"
#include "detray/io/common/homogeneous_material_writer.hpp"
#include "detray/io/common/io_interface.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/materials/material_slab.hpp"

// System include(s)
#include <string>

namespace detray {

/// @brief Abstract base class for material maps writers
template <class detector_t>
class material_map_writer
    : public detail::grid_writer<
          detector_t, material_slab<typename detector_t::scalar_type>> {

    using scalar_t = typename detector_t::scalar_type;
    using material_t = material_slab<scalar_t>;

    using base_type =
        detail::grid_writer<detector_t,
                            material_slab<typename detector_t::scalar_type>>;
    using grid_writer_t = base_type;
    using mat_writer_t = homogeneous_material_writer<detector_t>;

    protected:
    /// Tag the writer as "material_map"
    inline static const std::string tag = "material_maps";

    public:
    /// Same constructors for this class as for base_type
    using base_type::base_type;

    protected:
    /// Serialize the header information into its payload
    static auto write_header(const detector_t& det,
                             const std::string_view det_name) {

        return grid_writer_t::write_header(tag, det.material_store(), det_name);
    }

    /// Serialize the material description of a detector @param det into its io
    /// payload
    static detector_grids_payload<material_slab_payload> serialize(
        const detector_t& det, const typename detector_t::name_map&) {

        detector_grids_payload<material_slab_payload> grids_data;

        // How to serialize a material slab in the grid
        auto mat_serializer = [](const material_t& mat) {
            return mat_writer_t::serialize(mat, dindex_invalid);
        };

        for (const auto& vol_desc : det.volumes()) {

            /// Check if a surface has a metrial map
            for (const auto& sf_desc : det.surface_lookup()) {
                if (sf_desc.volume() != vol_desc.index()) {
                    continue;
                }

                const auto& mat_link = sf_desc.material();
                // Don't look at empty links
                if (mat_link.is_invalid() or
                    mat_link.id() == detector_t::materials::id::e_none) {
                    continue;
                }

                // Generate the payload
                grid_writer_t::serialize(det.material_store(), mat_link,
                                         sf_desc.index(), grids_data,
                                         mat_serializer);
            }
        }

        return grids_data;
    }
};

}  // namespace detray
