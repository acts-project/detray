/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/grid_axis.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/io/json/json.hpp"
#include "detray/io/json/json_algebra_io.hpp"
#include "detray/io/json/json_geometry_io.hpp"

// System include(s).
#include <array>
#include <optional>
#include <vector>

namespace detray {

void to_json(nlohmann::ordered_json& j, const axis_payload& a) {
    j["label"] = static_cast<unsigned int>(a.label);
    j["bounds"] = static_cast<unsigned int>(a.bounds);
    j["binning"] = static_cast<unsigned int>(a.binning);
    j["edges"] = a.edges;
    j["bins"] = a.bins;
}

void from_json(const nlohmann::ordered_json& j, axis_payload& a) {
    a.binning = static_cast<n_axis::binning>(j["binning"]);
    a.bounds = static_cast<n_axis::bounds>(j["bounds"]);
    a.label = static_cast<n_axis::label>(j["label"]);
    a.edges = j["edges"].get<std::vector<real_io>>();
    a.bins = j["bins"];
}

void to_json(nlohmann::ordered_json& j, const grid_payload& g) {
    nlohmann::ordered_json jaxes;
    for (const auto& a : g.axes) {
        jaxes.push_back(a);
    }
    j["axes"] = jaxes;
    j["entries"] = g.entries;
}

void from_json(const nlohmann::ordered_json& j, grid_payload& g) {
    nlohmann::ordered_json jaxes = j["axes"];
    for (auto jax : jaxes) {
        axis_payload a = jax;
        g.axes.push_back(a);
    }
    g.entries = j["entries"];
}
void to_json(nlohmann::ordered_json& j, const grid_objects_payload& g) {
    j["grid"] = g.grid;
    if (g.transform.has_value()) {
        j["transform"] = g.transform.value();
    }
}

void from_json(const nlohmann::ordered_json& j, grid_objects_payload& g) {
    g.grid = j["grid"];
    if (j.find("transform") != j.end()) {
        g.transform = j["transform"];
    }
}

void to_json(nlohmann::ordered_json& j, const links_payload& l) {
    nlohmann::ordered_json js;
    for (const auto& so : l.single_links) {
        js.push_back(so);
    }
    j["single_links"] = js;
    if (l.grid_links.has_value()) {
        j["grid_links"] = l.grid_links.value();
    }
}

void from_json(const nlohmann::ordered_json& j, links_payload& l) {
    nlohmann::ordered_json jsl = j["single_links"];
    for (auto jl : jsl) {
        single_link_payload sl = jl;
        l.single_links.push_back(sl);
    }
    if (j.find("grid_links") != j.end()) {
        l.grid_links = j["grid_links"];
    }
}

void to_json(nlohmann::ordered_json& j, const detector_payload& d) {
    j["name"] = d.name;
    if (not d.volumes.empty()) {
        nlohmann::ordered_json jvolumes;
        for (const auto& v : d.volumes) {
            jvolumes.push_back(v);
        }
        j["volumes"] = jvolumes;
        j["volume_grid"] = d.volume_grid;
    }
}

void from_json(const nlohmann::ordered_json& j, detector_payload& d) {
    d.name = j["name"];
    if (j.find("volumes") != j.end()) {
        for (auto jvolume : j["volumes"]) {
            d.volumes.push_back(jvolume);
        }
        d.volume_grid = j["volume_grid"];
    }
}

}  // namespace detray