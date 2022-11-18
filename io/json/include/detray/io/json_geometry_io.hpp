/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <array>

#include "detray/io/io_payload.hpp"
#include "detray/io/json_algebra_io.hpp"
#include "detray/io/json_defs.hpp"
#include "detray/io/json_material_io.hpp"
#include "detray/io/json_utilities_io.hpp"

/// @brief  The detray JSON I/O is written in such a way that it
/// can read/write ACTS files that are written with the Detray
/// JSON I/O extension
namespace detray {

void to_json(nlohmann::json& j, const mask_payload& m) {
    j["type"] = static_cast<unsigned int>(m.type);
    j["boundaries"] = m.boundaries;
}

void from_json(const nlohmann::json& j, mask_payload& m) {
    m.type = static_cast<mask_payload::mask_type>(j["type"]);
    m.boundaries = j["boundaries"].get<std::vector<real_io>>();
}

void to_json(nlohmann::json& j, const surface_payload& s) {
    j["transform"] = s.transform;
    j["mask"] = s.mask;
    j["geoID"] = s.gid;
    j["material"] = s.material;
}

void from_json(const nlohmann::json& j, surface_payload& s) {
    s.transform = j["transform"];
    s.mask = j["mask"];
    s.gid = j["geoID"];
    s.material = j["material"];
}

void to_json(nlohmann::json& j, const portal_payload& p) {
    j["surface"] = p.surface;
    j["volume_links"] = p.volume_links;
}

void from_json(const nlohmann::json& j, portal_payload& p) {
    p.surface = j["surface"];
    p.volume_links = j["volume_links"];
}

void to_json(nlohmann::json& j, const volume_bounds_payload& vb) {
    j["values"] = vb.values;
    j["type"] = static_cast<unsigned int>(vb.type);
}

void from_json(const nlohmann::json& j, volume_bounds_payload& vb) {
    vb.values = j["values"].get<std::vector<real_io>>();
    vb.type = static_cast<volume_bounds_payload::volume_bounds_type>(j["type"]);
}

void to_json(nlohmann::json& j, const volume_payload& v) {
    j["name"] = v.name;
    j["transform"] = v.transform;
    j["volume_bounds"] = v.volume_bounds;
    nlohmann::json pjson;
    for (const auto& p : v.portals) {
        pjson.push_back(p);
    }
    j["portals"] = pjson;
    if (not v.surfaces.empty()) {
        nlohmann::json sjson;
        for (const auto& s : v.surfaces) {
            sjson.push_back(s);
        }
        j["surfaces"] = sjson;
        j["surface_links"] = v.surface_links;
    }
}

void from_json(const nlohmann::json& j, volume_payload& v) {
    v.name = j["name"];
    v.transform = j["transform"];
    v.volume_bounds = j["volume_bounds"];
    for (auto jp : j["portals"]) {
        portal_payload p = jp;
        v.portals.push_back(p);
    }
    if (j.find("surfaces") != j.end()) {
        for (auto js : j["surfaces"]) {
            surface_payload s = js;
            v.surfaces.push_back(s);
        }
        v.surface_links = j["surface_links"];
    }
}

void to_json(nlohmann::json& j, const detector_payload& d) {
    j["name"] = d.name;
    if (not d.volumes.empty()) {
        nlohmann::json jvolumes;
        for (const auto& v : d.volumes) {
            jvolumes.push_back(v);
        }
        j["volumes"] = jvolumes;
        j["volume_grid"] = d.volume_grid;
    }
}

void from_json(const nlohmann::json& j, detector_payload& d) {
    d.name = j["name"];
    if (j.find("volumes") != j.end()) {
        for (auto jvolume : j["volumes"]) {
            d.volumes.push_back(jvolume);
        }
        d.volume_grid = j["volume_grid"];
    }
}

}  // namespace detray
