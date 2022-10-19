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
        const auto& n_volumes = det.volumes().size();

        for (std::size_t i = 0; i < n_volumes; i++) {
        }
    }

    std::vector<scalar_type> radiation_lengths = {};
};

}  // namespace detray