/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/io/common/io_interface.hpp"

// System include(s)
#include <ios>
#include <iostream>
#include <string>

namespace detray {

/// @brief Reader that constructs a covfie field from file input
template <class detector_t>
class covfie_file_reader : public reader_interface<detector_t> {

    using base_type = reader_interface<detector_t>;

    protected:
    /// Tag the reader as magnetic field reader
    inline static const std::string tag = "bfield";

    public:
    /// Same constructors for this class as for base_type
    using base_type::base_type;

    /// Covfie constructs the field directly from the input stream in
    /// @param file
    auto deserialize(io::detail::file_handle& file) const {
        return typename detector_t::bfield_type(*file);
    }
};

/// @brief Generic magnetic field reader
template <class detector_t, template <class> class reader_t>
class bfield_reader final : public reader_t<detector_t> {

    using base_reader = reader_t<detector_t>;

    public:
    /// Default file ending for covfie input
    bfield_reader() : base_reader{".cvf"} {}

    /// Read from file with given extension @param ext
    bfield_reader(const std::string& ext) : base_reader{ext} {}

    /// Reads the geometry from file with the given name
    virtual void read(detector_t& det, typename detector_t::name_map&,
                      const std::string& file_name) override {

        io::detail::file_handle file{file_name,
                                     std::ios_base::binary | std::ios_base::in};

        det.set_bfield(base_reader::deserialize(file));
    }
};

/// Read a covfie field from file and add it to a detray detector
template <typename detector_t>
using covfie_reader = bfield_reader<detector_t, covfie_file_reader>;

}  // namespace detray
