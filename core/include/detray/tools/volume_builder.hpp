/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/geometry.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/tools/volume_builder_interface.hpp"
// @TODO: Remove once material map writing becomes available
#include "detray/surface_finders/accelerator_grid.hpp"

// System include(s)
#include <memory>
#include <string>

namespace detray {

namespace detail {

/// A functor to update the mask index in surface descriptors
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

    /// Parametrized Constructor
    ///
    /// @param id flags the type of volume geometry (e.g. cylindrical, cuboid)
    /// @param idx the index of the volume in the detector volume container
    volume_builder(const volume_id id, const dindex idx = dindex_invalid)
        : m_has_accel{false}, m_volume{id} {

        m_volume.set_index(idx);

        // The first acceleration data structure in every volume is a brute
        // force method that will at least contain the portals
        m_volume
            .template set_link<static_cast<typename volume_type::object_id>(0)>(
                detector_t::accel::id::e_default, 0);
    };

    /// Adds the @param name of the volume to a @param name_map
    template <typename name_map>
    DETRAY_HOST void add_name(const name_map& names, std::string&& name) {
        names.at(vol_index()) = std::move(name);
    }

    /// @returns the volume index in the detector volume container
    DETRAY_HOST
    auto vol_index() -> dindex override { return m_volume.index(); }

    /// Toggles whether sensitive surfaces are added to the brute force method
    DETRAY_HOST
    void has_accel(bool toggle) override { m_has_accel = toggle; }

    /// @returns whether sensitive surfaces are added to the brute force method
    DETRAY_HOST
    bool has_accel() const override { return m_has_accel; }

    /// Access to the volume under construction - const
    DETRAY_HOST
    auto operator()() const -> const
        typename detector_t::volume_type& override {
        return m_volume;
    }

    /// Access to the volume under construction - non-const
    DETRAY_HOST
    auto operator()() -> typename detector_t::volume_type& override {
        return m_volume;
    }

    /// Build the volume with internal surfaces and portals and add it to the
    /// detector instance @param det
    DETRAY_HOST
    auto build(detector_t& det, typename detector_t::geometry_context ctx = {})
        -> typename detector_t::volume_type* override {
        // Prepare volume data
        m_volume.set_index(static_cast<dindex>(det.volumes().size()));

        m_volume.set_transform(det.transform_store().size());
        det.transform_store().push_back(m_trf);

        // Add all data from the builder to the detector containers
        add_to_detector(ctx, det);

        // Reset after the data was added to the detector
        m_surfaces.clear();
        m_transforms.clear(ctx);
        m_masks.clear_all();

        // Pass to decorator builders
        return &(det.volumes().back());
    }

    /// Adds a placement transform @param trf for the volume
    DETRAY_HOST
    void add_volume_placement(
        const typename detector_t::transform3& trf = {}) override {
        m_trf = trf;
    }

    /// Constructs a placement transform with identity rotation and translation
    /// @param t for the volume
    DETRAY_HOST
    void add_volume_placement(const typename detector_t::point3& t) override {
        m_trf = typename detector_t::transform3{t};
    }

    /// Constructs a placement transform from axes @param x and @param z
    /// and the translation @param t for the volume
    DETRAY_HOST
    void add_volume_placement(const typename detector_t::point3& t,
                              const typename detector_t::vector3& x,
                              const typename detector_t::vector3& z) override {
        m_trf = typename detector_t::transform3{t, z, x, true};
    }

    /// Add data for (a) new surface(s) to the builder
    DETRAY_HOST
    void add_surfaces(
        std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
        typename detector_t::geometry_context ctx = {}) override {
        (*sf_factory)(m_volume, m_surfaces, m_transforms, m_masks, ctx);
    }

    protected:
    /// @returns Access to the surface descriptor data
    typename detector_t::surface_container& surfaces() override {
        return m_surfaces;
    }

    /// @returns Access to the surface/volume transform data
    typename detector_t::transform_container& transforms() override {
        return m_transforms;
    }

    /// @returns Access to the surface mask data
    typename detector_t::mask_container& masks() override { return m_masks; }

    /// Adds a new full set of volume components (e.g. transforms or masks)
    /// to the global detector data stores and updates all links.
    ///
    /// @param ctx is the geometry_context of the call
    /// @param det is the detector instance that the volume should be added to
    ///
    /// @note can throw an exception if input data is inconsistent
    template <geo_obj_ids surface_id = static_cast<geo_obj_ids>(0)>
    DETRAY_HOST auto add_to_detector(
        const typename detector_t::geometry_context ctx,
        detector_t& det) noexcept(false) -> void {

        // Append transforms
        const auto trf_offset = det.transform_store().size(ctx);
        det.append_transforms(std::move(m_transforms), ctx);

        // Update mask and transform index of surfaces and set the
        // correct index of the surface in container
        auto sf_offset{static_cast<dindex>(det.surface_lookup().size())};
        std::size_t n_portals{0u};
        for (auto& sf_desc : m_surfaces) {

            const auto sf = surface{det, sf_desc};

            sf.template visit_mask<detail::mask_index_update>(sf_desc);
            sf_desc.set_volume(m_volume.index());
            sf_desc.update_transform(trf_offset);
            sf_desc.set_index(sf_offset++);

            det.add_surface_to_lookup(sf_desc);

            if (sf_desc.is_portal()) {
                ++n_portals;
            }
        }

        // Add portals to brute force navigation method
        if (m_has_accel) {
            typename detector_t::surface_container portals{};
            portals.reserve(n_portals);

            std::copy_if(m_surfaces.begin(), m_surfaces.end(),
                         std::back_inserter(portals),
                         [](auto& sf_desc) { return !sf_desc.is_sensitive(); });

            det.append_portals(std::move(portals));
        } else {
            // No acceleration structure: Add all surfaces to brute force method
            det.append_portals(std::move(m_surfaces));
        }

        // Update the surface link in the volume. In the volume builder, all
        // surfaces are filled into the default brute_force accelerator.
        // For the other accelerators (grid etc.) there need to be dedicated
        // builders
        constexpr auto default_acc_id{detector_t::accel::id::e_default};
        m_volume.template set_link<surface_id>(
            default_acc_id,
            det.accelerator_store().template size<default_acc_id>() - 1u);

        // Append masks
        det.append_masks(std::move(m_masks));

        // Finally, add the volume descriptor
        det.volumes().push_back(m_volume);
    }

    /// Whether the volume will get an acceleration structure
    bool m_has_accel{false};

    /// Volume descriptor of the volume under construction
    typename detector_t::volume_type m_volume{};
    /// Placement of the volume under construction
    typename detector_t::transform3 m_trf{};

    /// Data of conatined surfaces
    /// @{
    typename detector_t::surface_container m_surfaces{};
    typename detector_t::transform_container m_transforms{};
    typename detector_t::mask_container m_masks{};
    /// @}
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

/// TODO: Remove once the material builder is used everywhere
/// A functor to update the material index in surface objects
struct material_index_update {

    template <typename group_t, typename index_t, typename surface_t>
    DETRAY_HOST inline void operator()(
        [[maybe_unused]] const group_t& group,
        [[maybe_unused]] const index_t& /*index*/,
        [[maybe_unused]] surface_t& sf) const {
        if constexpr (!detail::is_grid_v<typename group_t::value_type>) {
            sf.update_material(static_cast<dindex>(group.size()));
        }
    }
};

}  // namespace detail

}  // namespace detray
