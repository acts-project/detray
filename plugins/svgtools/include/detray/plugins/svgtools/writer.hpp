/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/detector_writer.hpp"

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <vector>

namespace detray::svgtools {

/// @brief Writes a collection of svgs objects to a single file.
template <typename container_t>
inline void write_svg(const std::string& path, const container_t& svgs) {
    actsvg::svg::file file;
    for (const actsvg::svg::object& obj : svgs) {
        file.add_object(obj);
    }
    detray::io::detail::file_handle stream{path, "",
                                           std::ios::out | std::ios::trunc};
    *stream << file;
}

/// @brief Writes an svg object to a file.
/// @note To avoid conflict, the ids of the svg objects must be unique.
inline void write_svg(const std::string& path, const actsvg::svg::object& svg) {
    write_svg(path, std::array{svg});
}

/// @brief Writes an svg objects to a file.
/// @note To avoid conflict, the ids of the svg objects must be unique.
inline void write_svg(const std::string& path,
                      const std::initializer_list<actsvg::svg::object>& svgs) {
    std::vector<actsvg::svg::object> arg = svgs;
    write_svg(path, arg);
}

}  // namespace detray::svgtools