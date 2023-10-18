/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/io_interface.hpp"
#include "detray/tools/detector_builder.hpp"
#include "detray/tools/volume_builder.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <ios>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace detray::detail {

/// @brief A writer for multiple detector components.
///
/// The class aggregates a number of different writers that are created by a
/// manager class and calls them once the detector data should be written to
/// file. The respective writers are responsible for file handling etc.
///
/// @note the writers are unordered and must not depend on a write order!
template <class detector_t>
class detector_component_writers final {

    using writer_ptr_t = std::unique_ptr<writer_interface<detector_t>>;

    public:
    /// Default constructor
    detector_component_writers() = default;

    /// Create a new writer of type @tparam writer_t
    template <template <typename...> class writer_t, typename... Args,
              std::enable_if_t<std::is_base_of_v<writer_interface<detector_t>,
                                                 writer_t<detector_t, Args...>>,
                               bool> = true>
    void add() {
        add(std::make_unique<writer_t<detector_t, Args...>>());
    }

    /// Attach an existing writer via @param w_ptr to the writers
    void add(writer_ptr_t&& w_ptr) { m_writers.push_back(std::move(w_ptr)); }

    /// Writes the full detector data of @param det to file by calling the
    /// writers, while using the name map @param names for the detector
    virtual void write(const detector_t& det,
                       const typename detector_t::name_map& names,
                       const std::ios_base::openmode mode) {
        // We have to at least write a geometry
        assert(m_writers.size() != 0u &&
               "No writers registered! Need at least a geometry writer");

        // Call the write method on all optional writers
        std::for_each(m_writers.begin(), m_writers.end(),
                      [&det, &names, mode](writer_ptr_t& writer) {
                          writer->write(det, names, mode);
                      });
    }

    private:
    /// The writers registered for the detector: geoemtry (mandatory!) plus
    /// e.g. material, grids...)
    std::vector<writer_ptr_t> m_writers;
};

/// @brief A reader for multiple detector components.
///
/// The class aggregates a number of different readers that are created by a
/// manager class and calls them once the detector data should be read in from
/// file.
///
/// @note This class is still WIP! E.g. can it be collapsed with the components
/// writers ?
template <class detector_t>
class detector_component_readers final {

    using reader_ptr_t = std::unique_ptr<reader_interface<detector_t>>;

    public:
    /// Default constructor
    detector_component_readers() = default;

    /// Create a new reader of type @tparam reader_t
    template <template <typename...> class reader_t, typename... Args,
              std::enable_if_t<std::is_base_of_v<reader_interface<detector_t>,
                                                 reader_t<detector_t, Args...>>,
                               bool> = true>
    void add(const std::string& file_name) {
        add(std::make_unique<reader_t<detector_t, Args...>>(), file_name);
    }

    /// Attach an existing reader via @param w_ptr to the readers
    void add(reader_ptr_t&& w_ptr, const std::string& file_name) {
        m_readers[file_name] = std::move(w_ptr);
    }

    /// Reads the full detector into @param det by calling the readers, while
    /// using the name map @param names for to write the volume names.
    virtual void read(detector_builder<typename detector_t::metadata,
                                       volume_builder>& det_builder,
                      typename detector_t::name_map& names) {

        // We have to at least read a geometry
        assert(m_readers.size() != 0u &&
               "No readers registered! Need at least a geometry reader");

        // Call the read method on all readers
        for (const auto& [file_name, reader] : m_readers) {
            reader->read(det_builder, names, file_name);
        }
    }

    private:
    /// The readers registered for the detector: geoemtry (mandatory!) plus
    /// e.g. material, grids...)
    std::map<std::string, reader_ptr_t> m_readers;
};

}  // namespace detray::detail
