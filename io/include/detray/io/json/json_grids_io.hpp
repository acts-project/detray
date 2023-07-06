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

void to_json(nlohmann::ordered_json& j, const grid_header_payload& h) {
    j["version"] = h.version;
    j["detector"] = h.detector;
    j["date"] = h.date;
    j["tag"] = h.tag;
    j["no. grids"] = h.n_grids;
}

void from_json(const nlohmann::ordered_json& j, grid_header_payload& h) {
    h.version = j["version"];
    h.detector = j["detector"];
    h.date = j["date"];
    h.tag = j["tag"];
    h.n_grids = j["no. grids"];
}

void to_json(nlohmann::ordered_json& j, const axis_payload& a) {
    j["label"] = static_cast<unsigned int>(a.label);
    j["bounds"] = static_cast<unsigned int>(a.bounds);
    j["binning"] = static_cast<unsigned int>(a.binning);
    j["bins"] = a.bins;
    j["edges"] = a.edges;
}

void from_json(const nlohmann::ordered_json& j, axis_payload& a) {
    a.binning = static_cast<n_axis::binning>(j["binning"]);
    a.bounds = static_cast<n_axis::bounds>(j["bounds"]);
    a.label = static_cast<n_axis::label>(j["label"]);
    a.bins = j["bins"];
    a.edges = j["edges"].get<std::vector<real_io>>();
}

void to_json(nlohmann::ordered_json& j, const grid_bin_payload& g) {
    j["loc_index"] = g.loc_index;
    j["content"] = g.content;
}

void from_json(const nlohmann::ordered_json& j, grid_bin_payload& g) {
    g.loc_index = j["loc_index"].get<std::vector<unsigned int>>();
    g.content = j["content"].get<std::vector<std::size_t>>();
}

void to_json(nlohmann::ordered_json& j, const grid_payload& g) {
    j["type"] = static_cast<unsigned int>(g.type);
    j["index"] = g.index;

    nlohmann::ordered_json jaxes;
    for (const auto& a : g.axes) {
        jaxes.push_back(a);
    }
    j["axes"] = jaxes;

    nlohmann::ordered_json jbins;
    for (const auto& bin : g.bins) {
        jbins.push_back(bin);
    }
    j["bins"] = jbins;
}

void from_json(const nlohmann::ordered_json& j, grid_payload& g) {
    g.type = static_cast<grid_payload::grid_type>(j["type"]);
    g.index = j["index"];

    nlohmann::ordered_json jaxes = j["axes"];
    for (auto jax : jaxes) {
        axis_payload a = jax;
        g.axes.push_back(std::move(a));
    }

    nlohmann::ordered_json jbins = j["bins"];
    for (auto jbin : jbins) {
        grid_bin_payload b = jbin;
        g.bins.push_back(std::move(b));
    }
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
    if (j.find("volumes") != j.end()) {
        for (auto jvolume : j["volumes"]) {
            d.volumes.push_back(jvolume);
        }
        d.volume_grid = j["volume_grid"];
    }
}

void to_json(nlohmann::ordered_json& j, const detector_grids_payload& d) {
    if (not d.grids.empty()) {
        nlohmann::ordered_json jgrids;
        for (const auto& gr : d.grids) {
            jgrids.push_back(gr);
        }
        j["grids"] = jgrids;
    }
}

void from_json(const nlohmann::ordered_json& j, detector_grids_payload& d) {
    if (j.find("grids") != j.end()) {
        for (auto jgrid : j["grids"]) {
            grid_payload grp = jgrid;
            d.grids.push_back(std::move(grp));
        }
    }
}

}  // namespace detray
