/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/materials/material.hpp"
#include "detray/utils/ratio.hpp"

// System include(s)
#include <utility>

namespace detray {

template <typename... material_types>
struct material_composition {

    static_assert(ratio_sum<typename material_types::ratio...>::is_sum_one,
                  "Sumation of ratios should be equal to 1");

    material_composition(material_types... materials) {}
};

}  // namespace detray