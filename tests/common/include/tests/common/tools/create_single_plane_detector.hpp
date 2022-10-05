/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "tests/common/tools/detector_metadata.hpp"

namespace detray {

template <detector_registry::default_detector::mask_ids mask_id,
          detector_registry::default_detector::material_ids material_id,
          template <typename, std::size_t> class array_t = std::array,
          template <typename...> class tuple_t = std::tuple,
          template <typename...> class vector_t = vecmem::vector,
          template <typename...> class jagged_vector_t = vecmem::jagged_vector>
struct single_plane_detector_creator {

    using registry_type = detector_registry::default_detector;
    using detector_type =
        detector<registry_type, array_t, tuple_t, vector_t, jagged_vector_t>;
    using surface_type = typename detector_type::surface_container::value_type;
    using edge_type = typename surface_type::edge_type;
    using mask_link_type = typename surface_type::mask_link;
    using material_link_type = typename surface_type::material_link;

    single_plane_detector_creator(vecmem::memory_resource& resource)
        : m_detector(resource),
          m_surfaces(&resource),
          m_masks(resource),
          m_materials(resource),
          m_transforms(resource) {}

    template <typename... Args>
    inline void set_transform(Args&&... args) {
        m_transforms.emplace_back(ctx, args...);
    }

    template <typename... Args>
    inline void set_mask(Args&&... args) {
        m_masks.template add_value<mask_id>(args...,
                                            edge_type{0, dindex_invalid});
    }

    template <typename... Args>
    inline void set_material(Args&&... args) {
        m_materials.template add_value<material_id>(args...);
    }

    detector_type build() {
        const auto trf_index = m_transforms.size(ctx);
        mask_link_type mask_link{mask_id, m_masks.template size<mask_id>()};
        material_link_type material_link{
            material_id, m_materials.template size<material_id>()};
        m_surfaces.emplace_back(trf_index, mask_link, material_link, 0,
                                dindex_invalid, false);

        m_detector.new_volume({0., 0., 0., 0., -M_PI, M_PI});
        typename detector_type::volume_type& vol =
            m_detector.volume_by_index(0);

        m_detector.add_objects_per_volume(ctx, vol, m_surfaces, m_masks,
                                          m_materials, m_transforms);

        return std::move(m_detector);
    }

    detector_type m_detector;
    typename detector_type::context ctx{};
    typename detector_type::surface_container m_surfaces;
    typename detector_type::mask_container m_masks;
    typename detector_type::material_container m_materials;
    typename detector_type::transform_container m_transforms;
};

}  // namespace detray