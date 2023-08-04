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
#include "detray/io/json/json_common_io.hpp"

// System include(s).
#include <array>
#include <optional>
#include <vector>

namespace detray {

inline void to_json(nlohmann::ordered_json& j, const grid_header_payload& h) {
    j["common"] = h.common;

    if (h.sub_header.has_value()) {
        const auto& grid_sub_header = h.sub_header.value();
        j["grid_count"] = grid_sub_header.n_grids;
    }
}

inline void from_json(const nlohmann::ordered_json& j, grid_header_payload& h) {
    h.common = j["common"];

    if (j.find("grid_count") != j.end()) {
        h.sub_header.emplace();
        auto& grid_sub_header = h.sub_header.value();
        grid_sub_header.n_grids = j["grid_count"];
    }
}

inline void to_json(nlohmann::ordered_json& j, const axis_payload& a) {
    j["label"] = static_cast<unsigned int>(a.label);
    j["bounds"] = static_cast<unsigned int>(a.bounds);
    j["binning"] = static_cast<unsigned int>(a.binning);
    j["bins"] = a.bins;
    j["edges"] = a.edges;
}

inline void from_json(const nlohmann::ordered_json& j, axis_payload& a) {
    a.binning = static_cast<n_axis::binning>(j["binning"]);
    a.bounds = static_cast<n_axis::bounds>(j["bounds"]);
    a.label = static_cast<n_axis::label>(j["label"]);
    a.bins = j["bins"];
    a.edges = j["edges"].get<std::vector<real_io>>();
}

inline void to_json(nlohmann::ordered_json& j, const grid_bin_payload& g) {
    j["loc_index"] = g.loc_index;
    j["content"] = g.content;
}

inline void from_json(const nlohmann::ordered_json& j, grid_bin_payload& g) {
    g.loc_index = j["loc_index"].get<std::vector<unsigned int>>();
    g.content = j["content"].get<std::vector<std::size_t>>();
}

inline void to_json(nlohmann::ordered_json& j, const grid_payload& g) {
    j["volume_link"] = g.volume_link;
    j["acc_link"] = g.acc_link;

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

    if (g.transform.has_value()) {
        j["transform"] = g.transform.value();
    }
}

inline void from_json(const nlohmann::ordered_json& j, grid_payload& g) {
    g.volume_link = j["volume_link"];
    g.acc_link = j["acc_link"];

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

    if (j.find("transform") != j.end()) {
        g.transform.emplace();
        g.transform = j["transform"];
    }
}

inline void to_json(nlohmann::ordered_json& j,
                    const detector_grids_payload& d) {
    if (not d.grids.empty()) {
        nlohmann::ordered_json jgrids;
        for (const auto& gr : d.grids) {
            jgrids.push_back(gr);
        }
        j["grids"] = jgrids;
    }
}

inline void from_json(const nlohmann::ordered_json& j,
                      detector_grids_payload& d) {
    if (j.find("grids") != j.end()) {
        for (auto jgrid : j["grids"]) {
            grid_payload grp = jgrid;
            d.grids.push_back(std::move(grp));
        }
    }
}

}  // namespace detray
