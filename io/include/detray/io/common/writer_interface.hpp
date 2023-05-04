/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <ios>
#include <string>

namespace detray {

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
    virtual void write(const detector_t&, const typename detector_t::name_map&,
                       const std::ios_base::openmode) = 0;

    protected:
    /// Extension that matches the file format of the respective writer
    std::string m_file_extension;
};

}  // namespace detray
