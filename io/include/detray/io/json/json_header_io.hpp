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

namespace detray {

void to_json(nlohmann::ordered_json& j, const common_header_payload& h) {
    j["version"] = h.version;
    j["detector"] = h.detector;
    j["date"] = h.date;
    j["tag"] = h.tag;
}

void from_json(const nlohmann::ordered_json& j, common_header_payload& h) {
    h.version = j["version"];
    h.detector = j["detector"];
    h.date = j["date"];
    h.tag = j["tag"];
}

void to_json(nlohmann::ordered_json& j, const header_payload<bool>& h) {
    j["common"] = h.common;
    // Do write the optional subheader here, but in the dedicated serializers
}

void from_json(const nlohmann::ordered_json& j, header_payload<bool>& h) {
    h.common = j["common"];
    // Do not look at the optional subheader here
}

}  // namespace detray
