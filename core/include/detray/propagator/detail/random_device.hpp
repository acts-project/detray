/** Detray plugins library, part of the ACTS project
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <random>

namespace detray {

namespace stepping {

struct default_random_device {
    /// Random generator
    std::random_device rd{};
    std::mt19937_64 generator{rd()};
    void set_seed(const uint_fast64_t sd) { generator.seed(sd); }
};

}  // namespace stepping

}  // namespace detray