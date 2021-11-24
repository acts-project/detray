/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <iterator>
#include <utility>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/object_registry.hpp"

namespace detray {

// Minimalistic geometry type for toy geometry.
template <typename volume_t, typename object_t,
          template <typename...> class vector_t = dvector>
struct toy_geometry {
    // typedefs
    using volume_type = volume_t;
    using portal = object_t;
    using surface = object_t;
    using portal_container = vector_t<portal>;
    using portal_links = typename object_t::edge_links;
    using object_registry_type = toy_object_registry;

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    toy_geometry(vecmem::memory_resource& resource)
        : _volumes(&resource), _objects(&resource) {}

    /** Constructor from toy_geometry_data
     **/
    template <typename toy_geometry_data_t,
              std::enable_if_t<
                  !std::is_same_v<toy_geometry<volume_t, object_t, vector_t>,
                                  toy_geometry_data_t> &&
                      !std::is_base_of_v<vecmem::memory_resource,
                                         toy_geometry_data_t>,
                  bool> = true>
    DETRAY_DEVICE toy_geometry(toy_geometry_data_t& geometry_data)
        : _volumes(geometry_data._volumes_data),
          _objects(geometry_data._objects_data) {}

    DETRAY_HOST
    void add_volumes(vector_t<volume_t>&& volumes) {
        _volumes = std::move(volumes);
    }

    DETRAY_HOST
    void add_objects(vector_t<object_t>&& objects) {
        _objects = std::move(objects);
    }

    // data containers
    vector_t<volume_t> _volumes;
    vector_t<object_t> _objects;
};

template <typename toy_geometry_t>
struct toy_geometry_data {
    using volume_type = typename toy_geometry_t::volume_type;
    using object_type = typename toy_geometry_t::portal;

    toy_geometry_data(toy_geometry_t& geometry)
        : _volumes_data(vecmem::get_data(geometry._volumes)),
          _objects_data(vecmem::get_data(geometry._objects)) {}

    vecmem::data::vector_view<volume_type> _volumes_data;
    vecmem::data::vector_view<object_type> _objects_data;
};

template <typename volume_t, typename object_t,
          template <typename...> class vector_t>
inline toy_geometry_data<toy_geometry<volume_t, object_t, vector_t> > get_data(
    toy_geometry<volume_t, object_t, vector_t>& geometry) {
    return geometry;
}

}  // namespace detray