/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
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

namespace detail {

/// A functor to update the mask index in surface objects
struct mask_index_update;

}  // namespace detail

/// @brief Provides basic functionality to build detector volumes
template <typename detector_t>
class volume_builder : public volume_builder_interface<detector_t> {

    public:
    using scalar_type = typename detector_t::scalar_type;
    template <typename T, std::size_t N>
    using array_type = typename detector_t::template array_type<T, N>;
    using volume_type = typename detector_t::volume_type;
    using geo_obj_ids = typename detector_t::geo_obj_ids;

    volume_builder() = default;

    /// Adds the @param name of the volume to a name map
    template <typename name_map>
    DETRAY_HOST void add_name(const name_map& names, std::string&& name) {
        names.at(get_vol_index()) = std::move(name);
    }

    /// @brief Adds an array of @param bounds to a volume.
    DETRAY_HOST
    void init_vol(detector_t& det, const volume_id id) override {
        det.volumes().emplace_back(id);
        m_volume = &(det.volumes().back());
        m_volume->set_index(static_cast<dindex>(det.volumes().size()) - 1);
        m_volume->template set_link<
            static_cast<typename volume_type::object_id>(0)>(
            detector_t::sf_finders::id::e_default, 0);
    };

    /// @brief Adds an array of @param bounds to a volume.
    DETRAY_HOST
    void init_vol(detector_t& det, const volume_id id,
                  const array_type<scalar_type, 6>& bounds) override {
        det.volumes().emplace_back(id, bounds);
        m_volume = &(det.volumes().back());
        m_volume->set_index(static_cast<dindex>(det.volumes().size()) - 1);
        m_volume->template set_link<
            static_cast<typename volume_type::object_id>(0)>(
            detector_t::sf_finders::id::e_default, 0);
    };

    DETRAY_HOST
    auto get_vol_index() -> dindex override { return m_volume->index(); }

    DETRAY_HOST
    auto operator()() -> const typename detector_t::volume_type& override {
        return *m_volume;
    }

    DETRAY_HOST
    auto build(detector_t& det, typename detector_t::geometry_context ctx = {})
        -> typename detector_t::volume_type* override {
        add_objects_per_volume(ctx, det);
        m_surfaces.clear();
        m_transforms.clear(ctx);
        m_masks.clear_all();

        // Pass to decorator builders
        return m_volume;
    }

    DETRAY_HOST
    void add_portals(
        std::shared_ptr<surface_factory_interface<detector_t>> pt_factory,
        typename detector_t::geometry_context ctx = {}) override {
        (*pt_factory)(*m_volume, m_surfaces, m_transforms, m_masks, ctx);
    }

    DETRAY_HOST
    void add_sensitives(
        std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
        typename detector_t::geometry_context ctx = {}) override {
        (*sf_factory)(*m_volume, m_surfaces, m_transforms, m_masks, ctx);
    }

    DETRAY_HOST
    void add_passives(
        std::shared_ptr<surface_factory_interface<detector_t>> ps_factory,
        typename detector_t::geometry_context ctx = {}) override {
        (*ps_factory)(*m_volume, m_surfaces, m_transforms, m_masks, ctx);
    }

    protected:
    /// Add a new full set of detector components (e.g. transforms or volumes)
    /// according to given geometry_context.
    ///
    /// @param ctx is the geometry_context of the call
    /// @param det is the detector instance that the volume should be added to
    ///
    /// @note can throw an exception if input data is inconsistent
    template <geo_obj_ids surface_id = static_cast<geo_obj_ids>(0)>
    DETRAY_HOST auto add_objects_per_volume(
        const typename detector_t::geometry_context ctx,
        detector_t& det) noexcept(false) -> void {

        // Append transforms
        const auto trf_offset = det.transform_store().size(ctx);
        det.append_transforms(std::move(m_transforms), ctx);

        // Update mask and transform index of surfaces and set a
        // unique barcode (index of surface in container)
        auto sf_offset{static_cast<dindex>(det.surfaces().size())};
        for (auto& sf : m_surfaces) {
            det.mask_store().template visit<detail::mask_index_update>(
                sf.mask(), sf);
            sf.update_transform(trf_offset);
            sf.set_index(sf_offset++);
        }

        // Append surfaces
        det.append_surfaces(std::move(m_surfaces));

        // Update the surface link in the volume. In the volume builder, all
        // surfaces are filled into the default brute_force accelerator.
        // For the other accelerators (grid etc.) there need to be dedicated
        // builders
        constexpr auto default_acc_id{
            detector_t::sf_finders::id::e_brute_force};
        m_volume->template set_link<surface_id>(
            default_acc_id,
            det.surface_store().template size<default_acc_id>() - 1u);

        // Append masks
        det.append_masks(std::move(m_masks));
    }

    typename detector_t::volume_type* m_volume{};

    typename detector_t::surface_container_t m_surfaces{};
    typename detector_t::transform_container m_transforms{};
    typename detector_t::mask_container m_masks{};
};

namespace detail {

/// A functor to update the mask index in surface objects
struct mask_index_update {

    template <typename group_t, typename index_t, typename surface_t>
    DETRAY_HOST inline void operator()(const group_t& group,
                                       const index_t& /*index*/,
                                       surface_t& sf) const {
        sf.update_mask(static_cast<dindex>(group.size()));
    }
};

}  // namespace detail

}  // namespace detray
