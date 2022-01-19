/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

namespace detray {

// detector_registy for the hard-coded values
struct detector_registry {
    enum class default_detector { n_grids = 1 };
    enum class toy_detector { n_grids = 10 };
    enum class tml_detector { n_grids = 192 };
};

}  // namespace detray