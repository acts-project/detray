/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <array>
#include <optional>
#include <vector>

#include "detray/io/json_defs.hpp"

namespace detray {

/// @brief axis definition
struct axis_payload {
    /// axis lookup type
    enum class axis_lookup : unsigned int {
        x = 0u,
        y = 1u,
        z = 2u,
        r = 3u,
        phi = 4u
    };
    /// How the axis is done
    enum class axis_type : unsigned int { equidistant = 0u, variable = 1u };
    /// How the axis is bound
    enum class axis_bracket : unsigned int { bound = 0u, closed = 1u };

    axis_type type = axis_type::equidistant;
    axis_bracket bracket = axis_bracket::bound;
    axis_lookup lookup = axis_lookup::r;

    std::vector<real_io> borders = {};
    std::size_t bins = 0u;
};

void to_json(nlohmann::json& j, const axis_payload& a) {
    j["type"] = static_cast<unsigned int>(a.type);
    j["bracket"] = static_cast<unsigned int>(a.bracket);
    j["lookup"] = static_cast<unsigned int>(a.lookup);
    j["borders"] = a.borders;
    j["bins"] = a.bins;
}

void from_json(const nlohmann::json& j, axis_payload& a) {
    a.type = static_cast<axis_payload::axis_type>(j["type"]);
    a.bracket = static_cast<axis_payload::axis_bracket>(j["bracket"]);
    a.lookup = static_cast<axis_payload::axis_lookup>(j["lookup"]);
    a.borders = j["borders"].get<std::vector<real_io>>();
    a.bins = j["bins"];
}

/// @brief axis definition
struct grid_payload {
    std::vector<axis_payload> axes = {};
    std::vector<std::vector<unsigned int>> entries = {};
};

void to_json(nlohmann::json& j, const grid_payload& g) {
    nlohmann::json jaxes;
    for (const auto& a : g.axes) {
        jaxes.push_back(a);
    }
    j["axes"] = jaxes;
    j["entries"] = g.entries;
}

void from_json(const nlohmann::json& j, grid_payload& g) {
    nlohmann::json jaxes = j["axes"];
    for (auto jax : jaxes) {
        axis_payload a = jax;
        g.axes.push_back(a);
    }
    g.entries = j["entries"];
}

/// @brief single object
struct single_object_payload {
    unsigned int link;
};

void to_json(nlohmann::json& j, const single_object_payload& so) {
    j = so.link;
}

void from_json(const nlohmann::json& j, single_object_payload& so) {
    so.link = j;
}

/// @brief grid objects
struct grid_objects_payload {
    grid_payload grid;
    std::optional<transform_payload> transform;
};

void to_json(nlohmann::json& j, const grid_objects_payload& g) {
    j["grid"] = g.grid;
    if (g.transform.has_value()) {
        j["transform"] = g.transform.value();
    }
}

void from_json(const nlohmann::json& j, grid_objects_payload& g) {
    g.grid = j["grid"];
    if (j.find("transform") != j.end()) {
        g.transform = j["transform"];
    }
}

/// @brief navigation definition
struct links_payload {
    std::vector<single_object_payload> single_links;
    std::optional<grid_objects_payload> grid_links;
};

void to_json(nlohmann::json& j, const links_payload& l) {
    nlohmann::json js;
    for (const auto& so : l.single_links) {
        js.push_back(so);
    }
    j["single_links"] = js;
    if (l.grid_links.has_value()) {
        j["grid_links"] = l.grid_links.value();
    }
}

void from_json(const nlohmann::json& j, links_payload& l) {
    nlohmann::json jsl = j["single_links"];
    for (auto jl : jsl) {
        single_object_payload sl = jl;
        l.single_links.push_back(sl);
    }
    if (j.find("grid_links") != j.end()) {
        l.grid_links = j["grid_links"];
    }
}

}  // namespace detray
