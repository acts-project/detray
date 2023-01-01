/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/tools/bin_association.hpp"
#include "detray/tools/grid_factory.hpp"
#include "detray/tools/volume_builder.hpp"
#include "detray/tools/volume_builder_interface.hpp"

namespace detray {

// Forward declare default type
namespace detail {

struct fill_by_pos;

}  // namespace detail

/// @brief Build a grid of a certain shape.
///
/// Decorator class to a volume builder that adds a grid as the volumes
/// geometry accelerator structure.
template <typename detector_t, typename grid_t,
          typename bin_filler_t = detail::fill_by_pos,
          typename grid_factory_t = grid_factory_type<grid_t>>
class grid_builder : public volume_decorator<detector_t> {

    public:
    using scalar_type = typename detector_t::scalar_type;

    // TODO: nullptr can lead to exceptions, remove in the full implementation
    DETRAY_HOST
    grid_builder(std::unique_ptr<volume_builder_interface<detector_t>>
                     vol_builder = nullptr)
        : volume_decorator<detector_t>(std::move(vol_builder)) {}

    /// Should the passive surfaces be added to the grid ?
    void add_passives(bool is_add_passive = true) {
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

    /// Fill grid using a bin filling strategy @tparam bin_filler_t .
    /// This can also be called without a volume builder
    template <typename volume_type, typename... Args>
    DETRAY_HOST void fill_grid(
        const detector_t &det, const volume_type &vol,
        const typename detector_t::geometry_context ctx = {},
        const bin_filler_t bin_filler = {}, Args &&... /*args*/) {
        bin_filler(m_grid, det, vol, ctx);
    }

    /// Fill grid using a bin filling strategy @tparam bin_filler_t .
    /// This can also be called without a volume builder
    DETRAY_HOST
    void fill_grid(
        std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
        typename detector_t::geometry_context ctx = {},
        const bin_filler_t bin_filler = {}) {
        (*sf_factory)(volume_decorator<detector_t>::get_vol_index(), m_surfaces,
                      m_transforms, m_masks, ctx);
        bin_filler(m_grid, m_surfaces, m_transforms, m_masks, ctx);
    }

    /// Overwrite, to add the sensitives to the grid, instead of the surface vec
    DETRAY_HOST
    void add_sensitives(
        std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
        typename detector_t::geometry_context ctx = {}) override {
        fill_grid(std::move(sf_factory), ctx);
    }

    /// Overwrite, to add the passives to the grid, instead of the surface vec
    DETRAY_HOST
    void add_passives(
        std::shared_ptr<surface_factory_interface<detector_t>> ps_factory,
        typename detector_t::geometry_context ctx = {}) override {
        if (m_add_passives) {
            fill_grid(std::move(ps_factory), ctx);
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

        // Add the surafes that were filled into the grid directly to the
        // detector and update their links
        const auto trf_offset{det.transform_store().size(ctx)};
        // TODO: make the grid directly iterable
        for (std::size_t gbin{0}; gbin < m_grid.nbins(); ++gbin) {
            for (auto &sf : m_grid.at(gbin)) {
                det.mask_store().template visit<detail::mask_index_update>(
                    sf.mask(), sf);
                sf.update_transform(trf_offset);
            }
        }
        // Add transforms and masks to detector
        det.append_masks(std::move(m_masks));
        det.append_transforms(std::move(m_transforms));

        // Add the grid to the detector and link it to its volume
        constexpr auto gid{detector_t::sf_finders::template get_id<grid_t>()};
        det.surface_store().template push_back<gid>(m_grid);
        vol_ptr->set_link(gid, det.surface_store().template size<gid>() - 1);

        return vol_ptr;
    }

    /// @returns access to the new grid
    DETRAY_HOST
    const auto &operator()() const { return m_grid; }

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
          typename algebra_t = __plugin::transform3<detray::scalar>,
          template <typename, typename> class... binning_ts>
using grid_builder_type = grid_builder<
    detector_t,
    typename grid_factory_t<value_t, serializer_t, populator_impl_t,
                            algebra_t>::
        template grid_type<coordinate_axes<
            typename grid_shape_t::template axes<e_bounds, binning_ts...>, true,
            host_container_types, algebra_t>>,
    grid_factory_t<value_t, serializer_t, populator_impl_t, algebra_t>>;

namespace detail {

/// Fill a grid surface finder by bin association.
///
/// @param grid the grid that should be added
/// @param det the detector from which to get the surface placements
/// @param vol the volume the surface finder should be added to
/// @param ctx the geometry context
struct bin_associator {

    template <typename detector_t, typename volume_type, typename grid_t>
    DETRAY_HOST auto operator()(
        grid_t &grid, detector_t &det, const volume_type &vol,
        const typename detector_t::geometry_context ctx = {}) const -> void {
        // Fill the volumes surfaces into the grid
        this->operator()(grid, det.surfaces(vol), det.mask_store(),
                         det.transform_store(), ctx);
    }

    template <typename grid_t, typename surface_container,
              typename mask_container, typename transform_container>
    DETRAY_HOST auto operator()(
        grid_t &grid, const surface_container &surfaces,
        const transform_container &transforms, const mask_container &masks,
        const typename transform_container::context_type ctx = {}) const
        -> void {
        // Fill the surfaces into the grid
        bin_association(ctx, surfaces, transforms, masks, grid, {0.1, 0.1},
                        false);
    }
};

/// Fill a surface grid using the surface translation.
///
/// @param grid the grid that should be added
/// @param det the detector from which to get the surface placements
/// @param vol the volume the surface finder should be added to
/// @param ctx the geometry context
struct fill_by_pos {

    template <typename detector_t, typename volume_type, typename grid_t>
    DETRAY_HOST auto operator()(
        grid_t &grid, const detector_t &det, const volume_type &vol,
        const typename detector_t::geometry_context ctx = {}) const -> void {
        this->operator()(grid, det.surfaces(vol), det.transform_store(),
                         det.mask_store(), ctx);
    }

    template <typename grid_t, typename surface_container,
              typename mask_container, typename transform_container>
    DETRAY_HOST auto operator()(
        grid_t &grid, const surface_container &surfaces,
        const transform_container &transforms, const mask_container & /*masks*/,
        const typename transform_container::context_type /*ctx*/ = {}) const
        -> void {

        // Fill the volumes surfaces into the grid
        for (const auto &sf : surfaces) {
            if (sf.is_portal()) {
                continue;
            }
            const auto &sf_trf = transforms[sf.transform()];
            const auto &t = sf_trf.translation();
            const auto loc_pos = grid.global_to_local(
                __plugin::transform3<detray::scalar>{}, t, t);
            grid.populate(loc_pos, sf);
        }
    }
};

}  // namespace detail

}  // namespace detray