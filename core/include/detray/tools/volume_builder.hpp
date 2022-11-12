/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/surface.hpp"

// System include(s)
#include <string>

namespace detray {

/// @brief Provides functionality to build detector volumes
template <typename ID>
class volume_builder {

    // public:

    // virtual ~volume_builder() = default;

    /// Adds the @param name of the volume to a name map
    // void add_name(const name_map &names) = 0;

    /// @brief Adds an array of @param bounds to a volume.
    ///
    /// The bounds are interpreted differently for different volume shapes.
    // void add_bounds(const array_t<scalar_t, 6> &bounds) = 0;

    /// @brief Add a type of surface to the volume by ID.
    // template<ID surface_id, typename... Args>
    // dindex_range add_surface(Args&&... args) = 0;
};

/// @brief Decorator for the volume builder.
template <typename ID>
class volume_decorator : public volume_builder<ID> {

    /*public:
    volume_decorator() = delete;
    virtual ~volume_decorator() = default;

    volume_decorator(volume_builder<ID> * vol_builder) : m_builder{vol_builder}
    {}

    private:
    volume_builder<ID> *m_builder;*/
};

}  // namespace detray