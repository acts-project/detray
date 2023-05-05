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

void to_json(nlohmann::ordered_json& j,
             const homogeneous_material_header_payload& h) {
    j["version"] = h.version;
    j["detector"] = h.detector;
    j["date"] = h.date;
    j["tag"] = h.tag;
    j["no. slabs"] = h.n_slabs;
    j["no. rods"] = h.n_rods;
}

void from_json(const nlohmann::ordered_json& j,
               homogeneous_material_header_payload& h) {
    h.version = j["version"];
    h.detector = j["detector"];
    h.date = j["date"];
    h.tag = j["tag"];
    h.n_slabs = j["no. slabs"];
    h.n_rods = j["no. rods"];
}

void to_json(nlohmann::ordered_json& j, const material_payload& m) {
    j["params"] = m.params;
}

void from_json(const nlohmann::ordered_json& j, material_payload& m) {
    m.params = j["params"];
}

void to_json(nlohmann::ordered_json& j, const material_slab_payload& m) {
    j["type"] = material_slab_payload::material_type::slab;
    j["index"] = m.index;
    j["thickness"] = m.thickness;
    j["material"] = m.mat;
}

void from_json(const nlohmann::ordered_json& j, material_slab_payload& m) {
    m.type = material_slab_payload::material_type::slab;
    m.index = j["index"];
    m.thickness = j["thickness"];
    m.mat = j["material"];
}

void to_json(nlohmann::ordered_json& j,
             const detector_homogeneous_material_payload& d) {
    if (not d.mat_slabs.empty()) {
        nlohmann::ordered_json jmats;
        for (const auto& m : d.mat_slabs) {
            jmats.push_back(m);
        }
        j["material_slabs"] = jmats;
    }
    if (d.mat_rods.has_value() and not d.mat_rods->empty()) {
        nlohmann::ordered_json jmats;
        for (const auto& m : d.mat_rods.value()) {
            jmats.push_back(m);
        }
        j["material_rods"] = jmats;
    }
}

void from_json(const nlohmann::ordered_json& j,
               detector_homogeneous_material_payload& d) {
    if (j.find("material_slabs") != j.end()) {
        for (auto jmats : j["material_slabs"]) {
            material_slab_payload mslp = jmats;
            d.mat_slabs.push_back(mslp);
        }
    }
    if (j.find("material_rods") != j.end()) {
        d.mat_rods = {};
        for (auto jmats : j["material_rods"]) {
            material_slab_payload mslp = jmats;
            d.mat_rods->push_back(mslp);
        }
    }
}

}  // namespace detray
