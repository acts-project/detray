/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <iomanip>
#include <sstream>

namespace detray::detail {

/// Generate the filename for every distinct event
inline std::string get_event_filename(std::size_t event,
                                      const std::string& suffix) {
    std::stringstream stream;
    stream << "event";
    stream << std::setfill('0') << std::setw(9) << event;
    stream << suffix;
    return stream.str();
}

}  // namespace detray::detail