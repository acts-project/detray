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

    // TODO: nullptr can lead to exceptions, remove in the full implementation
    DETRAY_HOST
    grid_builder(
        std::unique_ptr<volume_builder_interface<detector_t>> vol_builder)
        : volume_decorator<detector_t>(std::move(vol_builder)) {}

    /// Should the passive surfaces be added to the grid ?
    void set_add_passives(bool is_add_passive = true) {
        m_add_passives = is_add_passive;
    }

    /// Delegate init call depending on @param span type
    template <typename grid_shape_t>
    DETRAY_HOST void init_grid(
        const mask<grid_shape_t> &bounds,
        const std::array<std::size_t, grid_t::Dim> &n_bins,
        const std::array<std::vector<scalar_type>, grid_t::Dim> &ax_bin_edges =
            {{}}) {
        init_impl(bounds, n_bins, ax_bin_edges,
                  typename grid_t::axes_type::bounds{},
                  typename grid_t::axes_type::binnings{});
    }

    /// Build the empty grid from axis parameters
    DETRAY_HOST void init_grid(
        const std::vector<scalar_type> &spans,
        const std::vector<std::size_t> &n_bins,
        const std::vector<std::vector<scalar_type>> &ax_bin_edges = {{}}) {

        m_grid = m_factory.template new_grid<typename grid_t::local_frame_type>(
            spans, n_bins, ax_bin_edges, typename grid_t::axes_type::bounds{},
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

    /// Add the volume and the grid to the detector @param det
    DETRAY_HOST
    auto build(detector_t &det, typename detector_t::geometry_context ctx = {})
        -> typename detector_t::volume_type * override {

        using surface_desc_t = typename detector_t::surface_type;

        // Add the surfaces (portals and/or passives) that are owned by the vol
        typename detector_t::volume_type *vol_ptr =
            volume_decorator<detector_t>::build(det, ctx);

        // Take the surfaces that should be filled into the grid out of the
        // brute force finder
        const auto vol_idx{vol_ptr->index()};
        constexpr auto bf_id{detector_t::sf_finders::id::e_brute_force};
        auto &bf_search = det.surface_store().template get<bf_id>();

        // Grid has not been filled previously, fill it automatically
        if (m_grid.size() == 0u) {
            std::vector<surface_desc_t> surfaces{};
            for (auto itr = bf_search.all().begin();
                 itr != bf_search.all().end();) {
                if (itr->volume() != vol_idx) {
                    continue;
                }
                if (itr->is_sensitive() or
                    (m_add_passives and itr->is_passive())) {
                    surfaces.push_back(*itr);
                    bf_search.erase(itr);
                } else {
                    ++itr;
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
            for (auto itr = bf_search.all().begin();
                 itr != bf_search.all().end();) {
                if (itr->volume() != vol_idx) {
                    ++itr;
                    continue;
                }
                if (itr->is_sensitive() or
                    (m_add_passives and itr->is_passive())) {
                    const auto vol = det.volume_by_index(itr->volume());
                    // The current volume is already built, so the surface
                    // interface is safe to use
                    const auto sf = surface{det, *itr};
                    const auto t = sf.center(ctx);
                    const auto loc_pos =
                        m_grid.global_to_local(vol.transform(), t, t);
                    auto bin_content = m_grid.search(loc_pos);

                    for (surface_desc_t &sf_desc : bin_content) {
                        // Find the correct surface and update all links
                        if (sf_desc.index() == sf.index()) {
                            sf_desc = *itr;
                        }
                    }

                    bf_search.erase(itr);
                } else {
                    ++itr;
                }
            }
        }

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

    grid_factory_t m_factory{};
    typename grid_t::template type<true> m_grid{};
    bin_filler_t m_bin_filler{};
    bool m_add_passives{false};
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
