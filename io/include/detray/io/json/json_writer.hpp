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
#include "detray/io/common/grid_writer.hpp"
#include "detray/io/common/homogeneous_material_writer.hpp"
#include "detray/io/json/json.hpp"
#include "detray/io/json/json_serializers.hpp"

// System include(s)
#include <cassert>
#include <ios>
#include <stdexcept>
#include <string>

namespace detray {

/// @brief Class that adds json functionality to common writer types.
///
/// Assemble the json writers from the common writer types, which serialize a
/// detector into the io payloads, and this class, which does the file
/// handling and provides the json stream. It also inlcudes the respective
/// @c to_json and @c from_json functions for the payloads ("json_serializers").
///
/// @note The resulting writer types will fulfill @c writer_interface through
/// the common writers they are being extended with
template <class detector_t, template <class> class common_writer_t>
class json_writer final : public common_writer_t<detector_t> {

    using base_writer = common_writer_t<detector_t>;

    public:
    /// File gets created with the json file extension
    json_writer() : base_writer(".json") {}

    /// Writes the geometry to file with a given name
    virtual std::string write(const detector_t &det,
                              const typename detector_t::name_map &names,
                              const std::ios_base::openmode mode) override {
        // Assert output stream
        assert(((mode == std::ios_base::out) or
                (mode == (std::ios_base::out | std::ios_base::trunc))) &&
               "Illegal file mode for json writer");

        // By convention the name of the detector is the first element
        std::string det_name = "";
        if (not names.empty()) {
            det_name = names.at(0);
        }

        // Create a new file
        std::string file_stem{det_name + "_" + base_writer::tag};
        io::detail::file_handle file{file_stem, this->m_file_extension, mode};

        // Write some general information
        nlohmann::ordered_json out_json;
        out_json["header"] = base_writer::write_header(det, det_name);

        // Write the detector data into the json stream by using the
        // conversion functions defined in "detray/io/json/json_io.hpp"
        out_json["data"] = base_writer::serialize(det, names);

        // Write to file
        *file << std::setw(4) << out_json << std::endl;

        return file_stem + this->m_file_extension;
    }
};

/// Write the tracking geometry to file in json format
template <typename detector_t>
using json_geometry_writer = json_writer<detector_t, geometry_writer>;

/// Write a simple material description to file in json format
template <typename detector_t>
using json_homogeneous_material_writer =
    json_writer<detector_t, homogeneous_material_writer>;

/// Write the detector grid collections to file in json format
template <typename detector_t>
using json_grid_writer = json_writer<detector_t, grid_writer>;

}  // namespace detray
