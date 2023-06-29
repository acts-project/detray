/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
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

    // TODO: nullptr can lead to exceptions, remove in the full implementation
    DETRAY_HOST
    grid_builder(std::unique_ptr<volume_builder_interface<detector_t>>
                     vol_builder = nullptr)
        : volume_decorator<detector_t>(std::move(vol_builder)) {}

    /// Should the passive surfaces be added to the grid ?
    void set_add_passives(bool is_add_passive = true) {
        m_add_passives = is_add_passive;
    }

    /// Delegate init call depending on @param span type
    template <typename axis_span_t>
    DETRAY_HOST void init_grid(
        const axis_span_t &span,
        const std::array<std::size_t, grid_t::Dim> &n_bins,
        const std::array<std::vector<scalar_type>, grid_t::Dim> &ax_bin_edges =
            {{}}) {
        init_impl(span, n_bins, ax_bin_edges,
                  typename grid_t::axes_type::bounds{},
                  typename grid_t::axes_type::binnings{});
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

    /// Overwrite, to add the sensitives to the grid, instead of the surface vec
    DETRAY_HOST
    void add_sensitives(
        std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
        typename detector_t::geometry_context ctx = {}) override {
        (*sf_factory)(volume_decorator<detector_t>::operator()(), m_surfaces,
                      m_transforms, m_masks, ctx);
    }

    /// Overwrite, to add the passives to the grid, instead of the surface vec
    DETRAY_HOST
    void add_passives(
        std::shared_ptr<surface_factory_interface<detector_t>> ps_factory,
        typename detector_t::geometry_context ctx = {}) override {
        if (m_add_passives) {
            (*ps_factory)(volume_decorator<detector_t>::operator()(),
                          m_surfaces, m_transforms, m_masks, ctx);
        } else {
            volume_decorator<detector_t>::add_passives(std::move(ps_factory),
                                                       ctx);
        }
    }

    /// Add the volume and the grid to the detector @param det
    DETRAY_HOST
    auto build(detector_t &det, typename detector_t::geometry_context ctx = {})
        -> typename detector_t::volume_type * override {
        // Add the surfaces (portals and/or passives) that are owned by the vol
        typename detector_t::volume_type *vol_ptr =
            volume_decorator<detector_t>::build(det, ctx);
        m_bin_filler(m_grid, detector_volume{det, *vol_ptr}, m_surfaces,
                     m_transforms, m_masks, ctx);

        // Add the surfaces that were filled into the grid directly to the
        // detector and update their links
        const auto trf_offset{det.transform_store().size(ctx)};
        for (auto &sf_desc : m_grid.all()) {
            const auto sf = surface{det, sf_desc};
            sf.template visit_mask<detail::mask_index_update>(sf_desc);
            sf_desc.update_transform(trf_offset);
        }

        // Add transforms and masks to detector
        det.append_masks(std::move(m_masks));
        det.append_transforms(std::move(m_transforms));

        // Add the grid to the detector and link it to its volume
        constexpr auto gid{detector_t::sf_finders::template get_id<grid_t>()};
        det.surface_store().template push_back<gid>(m_grid);
        vol_ptr->template set_link<
            detector_t::volume_type::object_id::e_sensitive>(
            gid, det.surface_store().template size<gid>() - 1);

        return vol_ptr;
    }

    /// @returns access to the new grid
    DETRAY_HOST
    auto &get() { return m_grid; }

    protected:
    /// Build the empty grid from a mask instance
    template <typename grid_shape_t, typename... axis_bounds,
              typename... binning_ts>
    DETRAY_HOST void init_impl(
        const mask<grid_shape_t> &bounds,
        const std::array<std::size_t, grid_t::Dim> n_bins,
        const std::array<std::vector<scalar_type>, grid_t::Dim> &ax_bin_edges,
        std::tuple<axis_bounds...>, std::tuple<binning_ts...>) {

        m_grid = m_factory.template new_grid<axis_bounds..., binning_ts...>(
            bounds, n_bins, ax_bin_edges);
    }

    /// Build the empty grid from axis parameters
    template <typename grid_shape_t, typename... axis_bounds,
              typename... binning_ts>
    DETRAY_HOST void init_impl(
        const std::vector<scalar_type> &spans,
        const std::vector<std::size_t> *n_bins,
        const std::vector<std::vector<scalar_type>> &ax_bin_edges,
        std::tuple<axis_bounds...>, std::tuple<binning_ts...>) {

        m_grid =
            m_factory
                .template new_grid<grid_shape_t, axis_bounds..., binning_ts...>(
                    spans, n_bins, ax_bin_edges);
    }

    grid_factory_t m_factory{};
    typename grid_t::template type<true> m_grid{};
    bin_filler_t m_bin_filler{};
    bool m_add_passives{false};

    // surfaces that are filled into the grid, but not the volume
    typename detector_t::surface_container_t m_surfaces{};
    typename detector_t::transform_container m_transforms{};
    typename detector_t::mask_container m_masks{};
};

/// Grid builder from single components
template <typename detector_t,
          template <typename, template <std::size_t> class, typename, typename>
          class grid_factory_t,
          typename grid_shape_t, typename value_t,
          template <std::size_t> class serializer_t, typename populator_impl_t,
          n_axis::bounds e_bounds = n_axis::bounds::e_closed,
          typename algebra_t = typename detector_t::transform3,
          template <typename, typename> class... binning_ts>
using grid_builder_type = grid_builder<
    detector_t,
    typename grid_factory_t<value_t, serializer_t, populator_impl_t,
                            algebra_t>::
        template grid_type<coordinate_axes<
            typename grid_shape_t::template axes<e_bounds, binning_ts...>, true,
            host_container_types, algebra_t>>,
    grid_factory_t<value_t, serializer_t, populator_impl_t, algebra_t>>;

}  // namespace detray
