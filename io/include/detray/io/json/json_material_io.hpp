/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/payloads.hpp"
#include "detray/io/json/json.hpp"

// System include(s)
#include <array>

/// @brief  The detray JSON I/O is written in such a way that it
/// can read/write ACTS files that are written with the Detray
/// JSON I/O extension
namespace detray {

void to_json(nlohmann::ordered_json& j, const material_slab_payload& m) {
    j["params"] = m.slab;
    j["density_eff"] = m.density_eff;
}

void from_json(const nlohmann::ordered_json& j, material_slab_payload& m) {
    m.slab = j["params"];
    m.density_eff = j["density_eff"];
}

}  // namespace detray
