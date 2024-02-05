/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/geometry/detector_volume.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/tools/bin_association.hpp"
#include "detray/tools/bin_fillers.hpp"
#include "detray/tools/grid_factory.hpp"
#include "detray/tools/surface_factory_interface.hpp"
#include "detray/tools/volume_builder.hpp"
#include "detray/tools/volume_builder_interface.hpp"

// System include(s)
#include <array>
#include <cassert>
#include <memory>
#include <vector>

namespace detray {

/// @brief Build a grid of a certain shape.
///
/// Decorator class to a volume builder that adds a grid as the volumes
/// geometry accelerator structure.
template <typename detector_t, typename grid_t,
          typename bin_filler_t = detail::fill_by_pos,
          typename grid_factory_t = grid_factory_type<grid_t>>
class grid_builder final : public volume_decorator<detector_t> {

    public:
    using scalar_type = typename detector_t::scalar_type;

    /// Use the grid builder stand-alone
    DETRAY_HOST
    grid_builder() : volume_decorator<detector_t>(nullptr) {}

    /// Decorate a volume with a grid
    DETRAY_HOST
    grid_builder(
        std::unique_ptr<volume_builder_interface<detector_t>> vol_builder)
        : volume_decorator<detector_t>(std::move(vol_builder)) {
        // The grid builder provides an acceleration structure to the
        // volume, so don't add sensitive surfaces to the brute force method
        if (this->m_builder) {
            this->m_builder->has_accel(true);
        }
    }

    /// Should the passive surfaces be added to the grid ?
    void set_add_passives(bool is_add_passive = true) {
        m_add_passives = is_add_passive;
    }

    /// Delegate init call depending on @param span type
    template <typename grid_shape_t>
    DETRAY_HOST void init_grid(
        const mask<grid_shape_t> &m,
        const std::array<std::size_t, grid_t::dim> &n_bins,
        const std::array<std::vector<scalar_type>, grid_t::dim> &ax_bin_edges =
            {{}}) {

        static_assert(
            std::is_same_v<typename grid_shape_t::template local_frame_type<
                               typename detector_t::transform3>,
                           typename grid_t::local_frame_type>,
            "Mask has incorrect shape");

        m_grid = m_factory.template new_grid<grid_t>(m, n_bins, ax_bin_edges);
    }

    /// Build the empty grid from axis parameters
    DETRAY_HOST void init_grid(
        const std::vector<scalar_type> &spans,
        const std::vector<std::size_t> &n_bins,
        const std::vector<std::vector<scalar_type>> &ax_bin_edges = {{}}) {

        m_grid =
            m_factory.template new_grid<grid_t>(spans, n_bins, ax_bin_edges);
    }

    /// Fill grid from existing volume using a bin filling strategy
    /// This can also be called without a volume builder
    template <typename volume_type, typename... Args>
    DETRAY_HOST void fill_grid(
        const detector_t &det, const volume_type &vol,
        const typename detector_t::geometry_context ctx = {},
        const bin_filler_t bin_filler = {}, Args &&... args) {

        bin_filler(m_grid, det, vol, ctx, args...);
    }

    /// Fill grid from externally provided surfaces - temporary solution until
    /// the volume builders can be deployed in the toy detector
    template <typename volume_type, typename surface_container_t,
              typename transform_container_t, typename mask_container_t,
              typename... Args>
    DETRAY_HOST void fill_grid(
        const volume_type &vol, const surface_container_t &surfaces,
        const transform_container_t &transforms, const mask_container_t &masks,
        const typename detector_t::geometry_context ctx = {},
        const bin_filler_t bin_filler = {}, Args &&... args) {

        bin_filler(m_grid, vol, surfaces, transforms, masks, ctx, args...);
    }

    /// Add the volume and the grid to the detector @param det
    DETRAY_HOST
    auto build(detector_t &det, typename detector_t::geometry_context ctx = {})
        -> typename detector_t::volume_type * override {

        using surface_desc_t = typename detector_t::surface_type;

        // Add the surfaces (portals and/or passives) that are owned by the vol
        typename detector_t::volume_type *vol_ptr =
            volume_decorator<detector_t>::build(det, ctx);

        // Find the surfaces that should be filled into the grid
        const auto vol = detector_volume{det, vol_ptr->index()};

        // Grid has not been filled previously, fill it automatically
        if (m_grid.size() == 0u) {

            std::vector<surface_desc_t> surfaces{};
            for (auto &sf_desc : vol.surfaces()) {

                if (sf_desc.is_sensitive() or
                    (m_add_passives and sf_desc.is_passive())) {
                    surfaces.push_back(sf_desc);
                }
            }

            this->fill_grid(
                detector_volume{det,
                                volume_decorator<detector_t>::operator()()},
                surfaces, det.transform_store(), det.mask_store(), ctx);
        } else {
            // The grid is prefilled with surface descriptors that contain the
            // correct surface indices per bin (e.g. from file IO).
            // Now add the rest of the linking information, which is only
            // available after the volume builder ran
            for (auto &sf_desc : vol.surfaces()) {

                if (sf_desc.is_sensitive() or
                    (m_add_passives and sf_desc.is_passive())) {
                    // The current volume is already built, so the surface
                    // interface is safe to use
                    const auto sf = surface{det, sf_desc};
                    const auto &sf_trf = sf.transform(ctx);
                    const auto t = sf_trf.point_to_global(sf.centroid());
                    const auto loc_pos = m_grid.project(vol.transform(), t, t);
                    auto &bin_content = m_grid.search(loc_pos);

                    for (surface_desc_t &sf_in_grid : bin_content) {
                        // Find the correct surface and update all links
                        if (sf_in_grid.index() == sf.index()) {
                            sf_in_grid = sf_desc;
                        }
                    }
                }
            }
        }

        // Add the grid to the detector and link it to its volume
        constexpr auto gid{detector_t::accel::template get_id<grid_t>()};
        det.accelerator_store().template push_back<gid>(m_grid);
        vol_ptr->template set_accel_link<
            detector_t::volume_type::object_id::e_sensitive>(
            gid, det.accelerator_store().template size<gid>() - 1);

        return vol_ptr;
    }

    /// @returns access to the new grid
    DETRAY_HOST
    auto &get() { return m_grid; }

    protected:
    grid_factory_t m_factory{};
    typename grid_t::template type<true> m_grid{};
    bin_filler_t m_bin_filler{};
    bool m_add_passives{false};
};

/// Grid builder from single components
template <typename detector_t,
          template <class, template <std::size_t> class, typename>
          class grid_factory_t,
          typename grid_shape_t, typename bin_t,
          template <std::size_t> class serializer_t,
          axis::bounds e_bounds = axis::bounds::e_closed,
          typename algebra_t = typename detector_t::transform3,
          template <typename, typename> class... binning_ts>
using grid_builder_type = grid_builder<
    detector_t,
    typename grid_factory_t<bin_t, serializer_t, algebra_t>::template grid_type<
        axes<grid_shape_t, e_bounds, binning_ts...>>,
    grid_factory_t<bin_t, serializer_t, algebra_t>>;

}  // namespace detray
