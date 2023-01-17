/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/core/detector.hpp"
#include "detray/core/detector_metadata.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/tools/volume_builder.hpp"
#include "detray/tools/volume_builder_interface.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// Covfie include(s)
#include <covfie/core/field.hpp>

// System include(s)
#include <memory>
#include <vector>

namespace detray {

/// @brief Provides functionality to build a detray detector volume by volume
///
/// @tparam metadata the type definitions for the detector
/// @tparam volume_builder_t the basic volume builder to be used for the
///                          geometry data
/// @tparam volume_data_t the data structure that hold the volume builders
template <typename metadata = default_metadata,
          template <typename> class volume_builder_t = volume_builder,
          template <typename...> class volume_data_t = std::vector>
class detector_builder {
    public:
    using detector_type =
        detector<metadata, covfie::field, host_container_types>;

    /// Empty detector builder
    detector_builder() = default;

    /// Add a new volume builder that will build a volume of the shape given by
    /// @param id
    template <typename... Args>
    DETRAY_HOST auto new_volume(const volume_id id, Args&&... args)
        -> volume_builder_interface<detector_type>* {

        m_volumes.push_back(std::make_unique<volume_builder_t<detector_type>>(
            id, static_cast<dindex>(m_volumes.size()),
            std::forward<Args>(args)...));

        return m_volumes.back().get();
    }

    /// Decorate a volume builder at position @param volume_idx with more
    /// functionality
    template <template <typename> class builder_t>
    DETRAY_HOST auto decorate(dindex volume_idx)
        -> volume_builder_interface<detector_type>* {

        m_volumes[volume_idx] = std::make_unique<builder_t<detector_type>>(
            std::move(m_volumes[volume_idx]));

        return m_volumes[volume_idx].get();
    }

    /// Access a particular volume builder by volume index @param volume_idx
    DETRAY_HOST
    auto operator[](dindex volume_idx)
        -> volume_builder_interface<detector_type>* {
        return m_volumes[volume_idx].get();
    }

    /// Assembles the final detector from the volumes builders and allocates
    /// the detector containers with the memory resource @param resource
    DETRAY_HOST
    auto build(vecmem::memory_resource& resource) -> detector_type {

        detector_type det{resource};

        for (auto& vol_builder : m_volumes) {
            vol_builder->build(det);
        }

        // TODO: Add sorting, data deduplication etc. here later...

        return det;
    }

    protected:
    /// Data structure that holds a volume builder for every detector volume
    volume_data_t<std::unique_ptr<volume_builder_interface<detector_type>>>
        m_volumes{};
};

}  // namespace detray
