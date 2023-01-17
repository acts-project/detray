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
#include "detray/io/json/json_algebra_io.hpp"
#include "detray/io/json/json_header_io.hpp"

// System include(s)
#include <array>
#include <optional>
#include <vector>

/// @brief  The detray JSON I/O is written in such a way that it
/// can read/write ACTS files that are written with the Detray
/// JSON I/O extension
namespace detray {

void to_json(nlohmann::ordered_json& j, const geo_header_payload& h) {
    j["common"] = h.common;

    if (h.sub_header.has_value()) {
        const auto& geo_sub_header = h.sub_header.value();
        j["volume_count"] = geo_sub_header.n_volumes;
        j["surface_count"] = geo_sub_header.n_surfaces;
    }
}

void from_json(const nlohmann::ordered_json& j, geo_header_payload& h) {
    h.common = j["common"];

    if (j.find("volume_count") != j.end() and
        j.find("surface_count") != j.end()) {
        h.sub_header.emplace();
        auto& geo_sub_header = h.sub_header.value();
        geo_sub_header.n_volumes = j["volume_count"];
        geo_sub_header.n_surfaces = j["surface_count"];
    }
}

void to_json(nlohmann::ordered_json& j, const single_link_payload& so) {
    j = so.link;
}

void from_json(const nlohmann::ordered_json& j, single_link_payload& so) {
    so.link = j;
}

void to_json(nlohmann::ordered_json& j, const mask_payload& m) {
    j["shape"] = static_cast<unsigned int>(m.shape);
    j["volume_link"] = m.volume_link;
    j["boundaries"] = m.boundaries;
}

void from_json(const nlohmann::ordered_json& j, mask_payload& m) {
    m.shape = static_cast<mask_payload::mask_shape>(j["shape"]);
    m.volume_link = j["volume_link"];
    m.boundaries = j["boundaries"].get<std::vector<real_io>>();
}

void to_json(nlohmann::ordered_json& j, const material_link_payload& m) {
    j["type"] = static_cast<unsigned int>(m.type);
    j["index"] = m.index;
}

void from_json(const nlohmann::ordered_json& j, material_link_payload& m) {
    m.type = static_cast<material_link_payload::material_type>(j["type"]);
    m.index = j["index"];
}

void to_json(nlohmann::ordered_json& j, const surface_payload& s) {
    j["barcode"] = s.barcode;
    j["type"] = static_cast<unsigned int>(s.type);
    j["source"] = s.source;
    j["transform"] = s.transform;
    j["mask"] = s.mask;
    if (s.material.has_value()) {
        j["material"] = s.material.value();
    }
}

void from_json(const nlohmann::ordered_json& j, surface_payload& s) {
    s.barcode = j["barcode"];
    s.type = static_cast<detray::surface_id>(j["type"]);
    s.source = j["source"];
    s.transform = j["transform"];
    s.mask = j["mask"];
    if (j.find("material") != j.end()) {
        s.material = j["material"];
    }
}

void to_json(nlohmann::ordered_json& j, const acc_links_payload& al) {
    j["type"] = static_cast<unsigned int>(al.type);
    j["index"] = al.index;
}

void from_json(const nlohmann::ordered_json& j, acc_links_payload& al) {
    al.type = static_cast<acc_links_payload::acc_type>(j["type"]);
    al.index = j["index"];
}

void to_json(nlohmann::ordered_json& j, const volume_payload& v) {
    j["name"] = v.name;
    j["index"] = v.index;
    j["type"] = v.type;
    j["transform"] = v.transform;
    nlohmann::ordered_json sjson;
    for (const auto& s : v.surfaces) {
        sjson.push_back(s);
    }
    j["surfaces"] = sjson;
    if (v.acc_links.has_value() and !v.acc_links.value().empty()) {
        nlohmann::ordered_json ljson;
        for (const auto& al : v.acc_links.value()) {
            ljson.push_back(al);
        }
        j["acc_links"] = ljson;
    }
}

void from_json(const nlohmann::ordered_json& j, volume_payload& v) {
    v.name = j["name"];
    v.index = j["index"];
    v.type = j["type"];
    v.transform = j["transform"];
    for (auto js : j["surfaces"]) {
        surface_payload s = js;
        v.surfaces.push_back(s);
    }
    if (j.find("acc_links") != j.end()) {
        v.acc_links = {};
        for (auto jl : j["acc_links"]) {
            acc_links_payload al = jl;
            v.acc_links->push_back(al);
        }
    }
}

}  // namespace detray
