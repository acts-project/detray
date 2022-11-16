/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <array>

#include "detray/io/json_algebra_io.hpp"
#include "detray/io/json_defs.hpp"

/// @brief  The detray JSON I/O is written in such a way that it
/// can read/write ACTS files that are written with the Detray
/// JSON I/O extension
namespace detray {

/// @brief  A payload object for masks
struct material_slab_payload {
    std::array<real_io, 5u> slab;
};

void to_json(nlohmann::json& j, const material_slab_payload& m) {
    j = m.slab;
}

void from_json(const nlohmann::json& j, material_slab_payload& m) {
    m.slab = j.get<std::array<real_io, 5u> >();
}

}  // namespace detray