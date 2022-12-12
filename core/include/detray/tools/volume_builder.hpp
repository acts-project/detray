/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/geometry.hpp"
#include "detray/tools/volume_builder_interface.hpp"

// System include(s)
#include <memory>
#include <string>

namespace detray {

/// @brief Provides basic functionality to build detector volumes
template <typename detector_t>
class volume_builder : public volume_builder_interface<detector_t> {

    public:
    using scalar_type = typename detector_t::scalar_type;
    template <typename T, std::size_t N>
    using array_type = typename detector_t::template array_type<T, N>;
    using geo_obj_ids = typename detector_t::geo_obj_ids;

    volume_builder() = default;

    /// Adds the @param name of the volume to a name map
    template <typename name_map>
    DETRAY_HOST void add_name(const name_map& names, std::string&& name) {
        names.at(get_vol_index()) = std::move(name);
    }

    /// @brief Adds an array of @param bounds to a volume.
    DETRAY_HOST
    void init_vol(detector_t& det, const volume_id id,
                  const array_type<scalar_type, 6>& bounds) override {
        det.volumes().emplace_back(id, bounds);
        m_volume = &(det.volumes().back());
        m_volume->set_index(det.volumes().size() - 1);
        m_volume->set_sf_finder(detector_t::sf_finders::id::e_default, 0);
    };

    DETRAY_HOST
    virtual auto get_vol_index() -> dindex override {
        return m_volume->index();
    }

    DETRAY_HOST
    auto build(detector_t& det, typename detector_t::geometry_context ctx = {})
        -> typename detector_t::volume_type* override {
        det.add_objects_per_volume(ctx, *m_volume, m_surfaces, m_masks,
                                   m_transforms);
        return m_volume;
    }

    DETRAY_HOST
    void add_portals(
        std::shared_ptr<surface_factory_interface<detector_t>> pt_factory,
        typename detector_t::geometry_context ctx = {}) override {
        dindex_range pt_range{(*pt_factory)(m_volume->index(), m_surfaces,
                                            m_transforms, m_masks, ctx)};
        m_volume->template update_obj_link<geo_obj_ids::e_portal>(pt_range);
    }

    DETRAY_HOST
    void add_sensitives(
        std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
        typename detector_t::geometry_context ctx = {}) override {
        dindex_range sf_range{(*sf_factory)(m_volume->index(), m_surfaces,
                                            m_transforms, m_masks, ctx)};
        m_volume->template update_obj_link<geo_obj_ids::e_sensitive>(sf_range);
    }

    DETRAY_HOST
    void add_passives(
        std::shared_ptr<surface_factory_interface<detector_t>> ps_factory,
        typename detector_t::geometry_context ctx = {}) override {
        dindex_range ps_range{(*ps_factory)(m_volume->index(), m_surfaces,
                                            m_transforms, m_masks, ctx)};
        m_volume->template update_obj_link<geo_obj_ids::e_passive>(ps_range);
    }

    protected:
    typename detector_t::volume_type* m_volume{};

    typename detector_t::surface_container m_surfaces{};
    typename detector_t::transform_container m_transforms{};
    typename detector_t::mask_container m_masks{};
};

}  // namespace detray