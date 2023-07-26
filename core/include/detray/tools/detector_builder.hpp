/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/core/detector.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/tools/volume_builder_interface.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// Covfie include(s)
#include <covfie/core/field.hpp>

// System include(s)
#include <memory>

namespace detray {

/// @brief Provides functionality to build a detray detector volume by volume
template <typename metadata, template <typename> class volume_builder_t,
          template <typename...> class volume_data_t = std::vector>
class detector_builder {
    public:
    using detector_type =
        detector<metadata, covfie::field, host_container_types>;

    /// Empty detector builder
    detector_builder() = default;

    /// Add a new volume builder
    template <typename... Args>
    DETRAY_HOST auto new_volume(const volume_id id, Args&&... args)
        -> volume_builder_interface<detector_type>* {
        m_volumes.push_back(std::make_unique<volume_builder_t<detector_type>>(
            id, static_cast<dindex>(m_volumes.size()),
            std::forward<Args>(args)...));

        return m_volumes.back().get();
    }

    /// Decorate a volume builder with more functionality
    template <template <typename> class builder_t>
    DETRAY_HOST auto decorate(dindex volume_idx)
        -> volume_builder_interface<detector_type>* {
        m_volumes[volume_idx] = std::make_unique<builder_t<detector_type>>(
            std::move(m_volumes[volume_idx]));

        return m_volumes[volume_idx].get();
    }

    /// Access a particular volume builder
    DETRAY_HOST
    auto operator[](dindex volume_idx)
        -> volume_builder_interface<detector_type>* {
        return m_volumes[volume_idx].get();
    }

    /// Assembles the final detector from the volumes builders
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
