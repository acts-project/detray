/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/surface_finders/grid/axis.hpp"
#include "detray/surface_finders/grid/populator.hpp"
#include "detray/surface_finders/grid/serializer.hpp"

namespace detray {

namespace detail {

/// @brief Helper type to assemble an multi-axis from shapes tuple and binnings
template <bool is_owning, typename containers, typename local_frame, typename,
          typename>
struct multi_axis_assembler;

/// @brief Specialized struct to extract axis shapes from a tuple
template <bool is_owning, typename containers, typename local_frame,
          typename... axis_shapes, typename... binning_ts>
struct multi_axis_assembler<is_owning, containers, local_frame,
                            dtuple<axis_shapes...>, dtuple<binning_ts...>> {

    static_assert(sizeof...(axis_shapes) == sizeof...(binning_ts),
                  "Number of axis shapes for this mask and given binning types "
                  "don't match!");

    using type =
        n_axis::multi_axis<is_owning, local_frame,
                           n_axis::single_axis<axis_shapes, binning_ts>...>;
};

}  // namespace detail

/// Typedef for easier construction from mask shapes
template <typename shape_t, bool is_owning = true,
          typename containers = host_container_types,
          typename algebra_t = __plugin::transform3<detray::scalar>>
using coordinate_axes = typename detail::multi_axis_assembler<
    is_owning, containers,
    typename shape_t::template coordinate_type<algebra_t>,
    typename shape_t::axes::types,
    typename shape_t::axes::template binning<
        containers, typename algebra_t::scalar_type>>::type;

/// A two-dimensional grid for object storage
///
/// @tparam populator_t  is a prescription what to do when a bin gets
/// pupulated, it broadcasts also the value type
/// @tparam serialzier_t  type of the serializer to the storage represenations
/*template<typename value_t, typename serializer_t,
         typename populator_t, typename containers, typename algebra>
class grid_builder {

    public:

    static constexpr bool is_owning = true;
    using value_type = value_t;
    using populator_type = populator_t<value_type>;
    using serializer_type = serializer_t;

    grid(vecmem::memory_resource &resource) : m_resource(resource) {}

    /// Single grid
    template <typename mask_t>
    DETRAY_HOST void grid(bin_data_view_t &&bin_data, ) {
        using axes_t = coordinate_system<true, containers, algebra, mask_t,
binning_ts...>;

        axes_t axes = construct_axes<true, mask_t>();
        using grid_t = grid<value_type, serializer_type, axes_t, true,
containers, algebra>;

        grid_t grid(m_resource);

        for(auto&& bin : bin_data) {
            grid.populate(bin.pos(), bin.value());
        }
    }

    /// Fill/populate from volume
    template <typename detector_t, typename volume_t>
    DETRAY_HOST void grid(const volume_t &vol, const detector_t &det) {
        for (const auto [sf, sf_idx]: enumerate(vol, det.surfaces())) {
            ...
        }
    }

    /// Fill/populate from surface
    template <typename detector_t, typename surface_t,
              template <typename, typename> class... binning_ts>
    DETRAY_HOST void grid(const surface_t &sf, const detector_t &det) {
        using mask_t = typename detector_t::mask_defs ...
        using axes_t = coordinate_system<true, containers, algebra, mask_t,
binning_ts...>;
    }


    /// add owning grid to detector collection
    template<template<bool ownership> class grid_t, typename detector_t>
    DETRAY_HOST
    void add_to_collection(grid_t<ownership> &&g, detector_t &det) {
        auto gr_coll = detail::get<grid_t<not is_owning>>(det.grid_store());

        gr_coll.append(std::forward<grid_t<is_owning>>(g));
    }

    /// @returns access to the serialized data
    DETRAY_HOST
    serialized_storage &finalize() { return _serialized_data; }

    /// @returns reference to the serialized data
    DETRAY_HOST
    const serialized_storage &tmp_data_store() const {
        return _serialized_data;
    }

    private:
    storage_type m_serialized_data{};
    populator_type m_populator{};
    serializer_type m_serializer{};
};*/

}  // namespace detray