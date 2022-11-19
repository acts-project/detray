/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/io/io_payload.hpp"
#include "detray/io/json_defs.hpp"

/// @brief  The detray JSON I/O is written in such a way that it
/// can read/write ACTS files that are written with the Detray
/// JSON I/O extension
namespace detray {

void to_json(nlohmann::json& j, const transform_payload& t) {
    j["translation"] = t.tr;
    j["rotation"] = t.rot;
}

void from_json(const nlohmann::json& j, transform_payload& t) {
    t.tr = j["translation"];
    t.rot = j["rotation"];
}

}  // namespace detray
