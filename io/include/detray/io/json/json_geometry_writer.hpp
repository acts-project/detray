/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/geometry_writer.hpp"
#include "detray/io/json/json.hpp"
#include "detray/io/json/json_serializers.hpp"

// System include(s)
#include <ios>
#include <string>

namespace detray {

/// @brief Class that writes a tracking geometry to json file
template <class detector_t>
class json_geometry_writer final : public geometry_writer<detector_t> {

    using base_writer = geometry_writer<detector_t>;

    public:
    /// File gets created with a fixed @param extension
    json_geometry_writer() : geometry_writer<detector_t>("json") {}

    /// Writes the geometry to file with a given name
    virtual void write(const detector_t &det,
                       const std::string &name) override {
        // Create a new file
        io::detail::file_handle file{name + "_geometry", this->m_file_extension,
                                     std::ios_base::out};

        // Write the detector geometry into the json stream
        nlohmann::ordered_json out_json = base_writer::serialize(det);

        // Write to file
        *file << std::setw(4) << out_json << std::endl;
    }
};

}  // namespace detray
