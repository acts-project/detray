/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <filesystem>
#include <string>
#include <system_error>

namespace detray::io {

/// Check if a given file path exists and generate it if not
inline auto create_path(const std::string& outdir) {

    auto path = std::filesystem::path(outdir);

    if (not std::filesystem::exists(path)) {
        std::error_code err;
        if (!std::filesystem::create_directories(path, err)) {
            throw std::runtime_error(err.message());
        }
    }

    return path;
}

}  // namespace detray::io
