/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/detail/utils.hpp"
#include "detray/io/common/payloads.hpp"
#include "detray/tools/detector_builder.hpp"

// System include(s)
#include <ios>
#include <string>
#include <string_view>

namespace detray {

/// @brief Abstract base class for detray detector components readers
template <class detector_t>
class reader_interface {

    public:
    /// All readers must define a file extension
    reader_interface() = delete;

    /// Only accept files with a fixed @param extension
    reader_interface(const std::string& ext) : m_file_extension{ext} {}

    /// Default destructor
    virtual ~reader_interface() = default;

    /// Reads the respective detector component from file. Since the detector
    /// does not keep the volume names, the name map is also passed and
    /// filled.
    virtual void read(
        detector_builder<typename detector_t::metadata, volume_builder>&,
        typename detector_t::name_map&, const std::string&) = 0;

    protected:
    /// Extension that matches the file format of the respective reader
    std::string m_file_extension;
};

/// @brief Abstract base class for detray detector components writers
template <class detector_t>
class writer_interface {

    public:
    /// All writers must define a file extension
    writer_interface() = delete;

    /// File gets created with a fixed @param extension
    writer_interface(const std::string& ext) : m_file_extension{ext} {}

    /// Default destructor
    virtual ~writer_interface() = default;

    /// Writes the respective detector component to file. Since the detector
    /// does not provide the volume names, the name map is also passed.
    virtual std::string write(const detector_t&,
                              const typename detector_t::name_map&,
                              const std::ios_base::openmode) = 0;

    protected:
    /// Serialize the common header information using the detector name
    /// @param det_name and the file tag @param tag that describes the data file
    /// content
    static common_header_payload serialize(const std::string_view det_name,
                                           const std::string_view tag) {
        common_header_payload header_data;

        header_data.version = detail::get_detray_version();
        header_data.detector = det_name;
        header_data.tag = tag;
        header_data.date = detail::get_current_date();

        return header_data;
    }

    /// Extension that matches the file format of the respective writer
    std::string m_file_extension;
};

}  // namespace detray
