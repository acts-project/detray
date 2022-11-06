/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/masks/masks.hpp"
#include "detray/tools/grid_factory.hpp"
#include "detray/tools/volume_builder.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <iostream>
#include <string>

namespace detray {

namespace detail {
/// Fill a grid surface finder by bin association, then add it to the
/// detector.
///
/// New surface finder id can be been given explicitly. That is helpful, if
/// multiple sf finders have the same type in the tuple container. Otherwise
/// it is determined automatically.
///
/// @param ctx the geometry context
/// @param vol the volume the surface finder should be added to
/// @param grid the grid that should be added
struct bin_assoiator {

    template <typename detector_t, typename grid_t>
    auto operator()(grid_t &grid, detector_t &det,
                    const typename detector_t::volume_type &vol,
                    const typename detector_t::geometry_context ctx = {}) const
        -> void {
        // Fill the volumes surfaces into the grid
        bin_association(ctx, det, vol, grid, {0.1, 0.1}, false);
    }
};

/// Fill a grid surface finder by bin association, then add it to the
/// detector.
///
/// New surface finder id can be been given explicitly. That is helpful, if
/// multiple sf finders have the same type in the tuple container. Otherwise
/// it is determined automatically.
///
/// @param ctx the geometry context
/// @param vol the volume the surface finder should be added to
/// @param grid the grid that should be added
struct fill_by_pos {

    template <typename detector_t, typename grid_t>
    auto operator()(grid_t &grid, const detector_t &det,
                    const typename detector_t::volume_type &vol,
                    const typename detector_t::geometry_context /*ctx*/ = {})
        const -> void {

        // Fill the volumes surfaces into the grid
        const auto &trf_store = det.transform_store();
        for (const auto &[idx, sf] :
             detray::views::enumerate(det.surfaces(), vol)) {
            if (sf.is_portal()) {
                continue;
            }
            const auto &sf_trf = trf_store[sf.transform()];
            const auto &t = sf_trf.translation();
            const auto loc_pos = grid.global_to_local(
                __plugin::transform3<detray::scalar>{}, t, t);
            grid.populate(loc_pos, idx);
        }
    }
};

}  // namespace detail

/// @brief Build a grid of a certain shape.
///
/// Decorator class to a volume builder that adds a grid as the volumes
/// geometry accelerator structure.
template <typename ID, typename grid_t,
          typename grid_factory_t = grid_factory_type<grid_t>,
          typename scalar_t = scalar>
class grid_builder : public volume_decorator<ID> {

    public:
    grid_builder() = default;

    /// Delegate init call
    template <typename axis_span_t>
    void init(const axis_span_t &span,
              const std::array<std::size_t, grid_t::Dim> &n_bins,
              const std::array<std::vector<scalar_t>, grid_t::Dim>
                  &ax_bin_edges = {{}}) {
        /// m_volume_builder->init(args...);
        init_impl(span, n_bins, ax_bin_edges,
                  typename grid_t::axes_type::bounds{},
                  typename grid_t::axes_type::binnings{});
    }

    template <typename bin_filler_t, typename detector_t, typename... Args>
    void fill(const bin_filler_t &bin_filler, const detector_t &det,
              const typename detector_t::volume_type &vol,
              const typename detector_t::geometry_context ctx = {},
              Args &&.../*args*/) {
        /// m_volume_builder->fill(args...);
        bin_filler(m_grid, det, vol, ctx);
    }

    const auto &operator()() const { return m_grid; }

    private:
    template <typename grid_shape_t, typename... axis_bounds,
              typename... binning_ts>
    void init_impl(
        const mask<grid_shape_t> &bounds,
        const std::array<std::size_t, grid_t::Dim> n_bins,
        const std::array<std::vector<scalar_t>, grid_t::Dim> &ax_bin_edges,
        std::tuple<axis_bounds...>, std::tuple<binning_ts...>) {

        m_grid = m_factory.template new_grid<axis_bounds..., binning_ts...>(
            bounds, n_bins, ax_bin_edges);
    }

    template <typename grid_shape_t, typename... axis_bounds,
              typename... binning_ts>
    void init_impl(const std::vector<scalar_t> &spans,
                   const std::vector<std::size_t> *n_bins,
                   const std::vector<std::vector<scalar_t>> &ax_bin_edges,
                   std::tuple<axis_bounds...>, std::tuple<binning_ts...>) {

        m_grid =
            m_factory
                .template new_grid<grid_shape_t, axis_bounds..., binning_ts...>(
                    spans, n_bins, ax_bin_edges);
    }

    volume_decorator<ID> *m_volume_builder{nullptr};
    grid_factory_t m_factory{};
    typename grid_t::template type<true> m_grid{};
};

/// Grid builder from single components
template <typename ID,
          template <typename, template <std::size_t> class, typename, typename>
          class grid_factory_t,
          typename grid_shape_t, typename value_t,
          template <std::size_t> class serializer_t, typename populator_impl_t,
          n_axis::bounds e_bounds = n_axis::bounds::e_closed,
          typename algebra_t = __plugin::transform3<detray::scalar>,
          template <typename, typename> class... binning_ts>
using grid_builder_type = grid_builder<
    ID,
    typename grid_factory_t<value_t, serializer_t, populator_impl_t,
                            algebra_t>::
        template grid_type<coordinate_axes<
            typename grid_shape_t::template axes<e_bounds, binning_ts...>, true,
            host_container_types, algebra_t>>,
    grid_factory_t<value_t, serializer_t, populator_impl_t, algebra_t>,
    typename algebra_t::scalar_type>;

}  // namespace detray