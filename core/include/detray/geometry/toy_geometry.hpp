/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>
#include <utility>

#include "detray/definitions/detray_qualifiers.hpp"
#include "detray/geometry/object_registry.hpp"

namespace detray {

// Minimalistic geometry type for toy geometry.
template <typename volume_t, typename object_t,
          template <typename...> class vector_t = dvector>
struct toy_geometry {
    // typedefs
    using volume_type = volume_t;
    using portal = object_t;
    using portal_container = vector_t<portal>;
    using portal_links = typename object_t::edge_links;
    using object_registry_type = toy_object_registry;

    toy_geometry(vector_t<volume_t>&& volumes, vector_t<object_t>&& objects)
        : _volumes(volumes), _objects(objects) {}

    // data containers
    vector_t<volume_t> _volumes;
    vector_t<object_t> _objects;
};

}  // namespace detray