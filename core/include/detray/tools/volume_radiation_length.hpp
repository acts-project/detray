/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <climits>

namespace detray {

template <typename detector_t>
class volume_radiation_length {

    using scalar_type = detector_t::scalar_type;

    volume_radiation_length(const detector_t& det) {

        for (const auto& vol : det.volumes()) {
            const scalar_type vol_size = vol.volume_size();

            for (const auto [obj_idx, obj] : detray::views::enumerate(det->surfaces(), vol){



            }
        }
    }

    std::vector<scalar_type> radiation_lengths = {};
};

}  // namespace detray