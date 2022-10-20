/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/utils/ranges.hpp"

// System include(s).
#include <climits>

namespace detray {

template <typename detector_t>
class volume_radiation_length {

    using scalar_type = typename detector_t::scalar_type;

    struct get_mask_area {};

    struct get_material_area {};

    struct get_surface_rad_len {
        using output_type = scalar_type;

        template <typename material_group_t, typename index_t,
                  typename surface_t>
        DETRAY_HOST_DEVICE inline output_type operator()(
            const material_group_t &material_group, const index_t & /*index*/,
            const surface_t &surface) const {

            const auto &material_range = surface.material_range();

            for (const auto &mat :
                 detray::ranges::subrange(material_group, material_range)) {
            }

            return 0;
        }
    };

    volume_radiation_length(const detector_t &det) {

        // @todo: Consider the case where the volume is filled with gas

        const auto &mat_store = det.material_store();

        for (const auto &vol : det.volumes()) {
            const scalar_type vol_size = vol.volume_size();

            for (const auto [obj_idx, obj] :
                 detray::views::enumerate(det->surfaces(), vol)) {
            }
        }
    }

    std::vector<scalar_type> radiation_lengths = {};
};

}  // namespace detray