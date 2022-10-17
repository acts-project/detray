/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/accessor.hpp"  // detail::first_t
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/surface_finders/grid/detail/axis_binning.hpp"
#include "detray/surface_finders/grid/detail/axis_helpers.hpp"
#include "detray/surface_finders/grid/detail/axis_shape.hpp"
#include "detray/utils/type_registry.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <cstddef>
#include <type_traits>

namespace detray {

namespace n_axis {

/// @brief A single axis.
///
/// The axis ties together shape and binning behaviour with the bin edges
/// storage. The shape determines how bin indices are mapped at the over-
/// and underflow bins. The type of binning determines whether the axis has
/// regular or irregular binning. The bin edges needed to find bin indices
/// are not owned by the axis, but are passed to the binning type.
template <typename shape_t, typename binning_t>
struct single_axis {

    /// Make axis shape accessible
    using shape_type = shape_t;
    /// Make axis binning type accessible
    using binning_type = binning_t;

    // Extract container types
    using scalar_type = typename binning_type::scalar_type;
    using container_types = typename binning_type::container_types;
    template <typename T>
    using vector_type = typename binning_type::template vector_type<T>;
    template <typename T, std::size_t N>
    using array_type = typename binning_type::template array_type<T, N>;

    /// Defines the geometrical shape of the axis as a service:
    /// open, closed or circular
    shape_type m_shape{};
    /// Defines the binning on the axis as a service: regular vs irregular
    binning_type m_binning{};

    /// Parameterized constructor - empty binning
    single_axis() = default;

    /// Parametrized constructor that builds the binning scheme
    template <typename... Args>
    DETRAY_HOST_DEVICE single_axis(const dindex_range *indx_range,
                                   const vector_type<scalar_type> *edges)
        : m_shape(), m_binning(indx_range, edges) {}

    /// @returns the axis shape label, i.e. x, y, z, r or phi axis.
    DETRAY_HOST_DEVICE
    constexpr auto label() const -> n_axis::label { return shape_type::label; }

    /// @returns the type of shape of the axis, i.e. closed, open or circular.
    DETRAY_HOST_DEVICE
    constexpr auto shape() const -> n_axis::shape { return shape_type::type; }

    /// @returns the type of binning of the axis, i.e. regular or irregular.
    DETRAY_HOST_DEVICE
    constexpr auto binning() const -> n_axis::binning {
        return binning_type::type;
    }

    /// @returns the total number of bins
    DETRAY_HOST_DEVICE
    inline constexpr std::size_t nbins() const { return m_binning.nbins(); }

    /// @returns the width of a bin
    template <typename... Args>
    DETRAY_HOST_DEVICE inline constexpr scalar_type bin_width(
        Args &&...args) const {
        return m_binning.bin_width(std::forward<Args &&>(args)...);
    }

    /// Given a value on the axis, find the correct bin.
    ///
    /// @param v is the value for the bin search
    ///
    /// @returns the bin index.
    DETRAY_HOST_DEVICE
    inline dindex bin(const scalar_type v) const {
        return m_shape.map(m_binning.bin(v), m_binning.nbins());
    }

    /// Given a value on the axis and a neighborhood, find the correct bin range
    ///
    /// @param v is the value for the bin search
    /// @param nhood is the neighborhood range (in #bins or value interval)
    ///
    /// @returns a dindex_range around the bin index.
    template <typename neighbor_t>
    DETRAY_HOST_DEVICE dindex_range
    range(const scalar_type v, const array_type<neighbor_t, 2> &nhood) const {
        return m_shape.map(m_binning.range(v, nhood), m_binning.nbins());
    }

    /// @returns the bin edges for a given @param ibin .
    DETRAY_HOST_DEVICE
    inline array_type<scalar_type, 2> bin_edges(const dindex ibin) const {
        return m_binning.bin_edges(ibin);
    }

    /// @returns the values of the bin edges. Is a succession of lower edges.
    DETRAY_HOST_DEVICE
    inline vector_type<scalar_type> bin_edges() const {
        return m_binning.bin_edges();
    }

    /// @returns the axis span [min, max).
    DETRAY_HOST_DEVICE
    inline array_type<scalar_type, 2> span() const { return m_binning.span(); }

    /// @returns the axis span [min, max).
    DETRAY_HOST_DEVICE
    inline scalar_type min() const { return m_binning.span()[0]; }

    /// @returns the axis span [min, max).
    DETRAY_HOST_DEVICE
    inline scalar_type max() const { return m_binning.span()[1]; }
};

/// @brief A collection of single axes.
///
/// Given a point in the grid local coordinate system, which is spanned by the
/// axes in this multi-axes type, the corresponding bin multi-index or
/// multi-index range is returned.
///
/// @note can be owning the data (as member of a standalone grid) or can be
/// non-owning if the grid is part of a larger collection.
template <bool ownership, typename local_frame_t, typename... axis_ts>
class multi_axis {

    public:
    /// Dimension of the local coordinate system that is spanned by the axes
    static constexpr dindex Dim = sizeof...(axis_ts);
    static constexpr bool is_owning = ownership;

    // Extract container types
    using scalar_type =
        typename detray::detail::first_t<axis_ts...>::scalar_type;
    using container_types =
        typename detray::detail::first_t<axis_ts...>::container_types;
    template <typename T>
    using vector_type = typename container_types::template vector_type<T>;
    template <typename... T>
    using tuple_type = typename container_types::template tuple_type<T...>;

    /// Projection onto local coordinate system that is spanned by the axes
    using local_frame_type = local_frame_t;

    /// Axes boundary/bin edges storage
    using boundary_storage_type = vector_type<dindex_range>;
    using edges_storage_type = vector_type<scalar_type>;

    /// Vecmem based multi-axis view type
    using view_type =
        dmulti_view<dvector_view<dindex_range>, dvector_view<scalar_type>>;
    using const_view_type = dmulti_view<dvector_view<const dindex_range>,
                                        dvector_view<const scalar_type>>;
    /// Interal storage type depends on whether the class owns the data or not
    using storage_type = std::conditional_t<
        is_owning, detail::multi_axis_data<container_types, scalar_type>,
        detail::multi_axis_view<container_types, scalar_type>>;

    /// Match an axis to its label at compile time
    using axis_reg = type_registry<n_axis::label, axis_ts...>;
    template <n_axis::label L>
    using label_matcher = typename axis_reg::template get_type<L, tuple_type>;

    /// Default constructor
    multi_axis() = default;

    /// Constructor with specific vecmem memory resource if the class owns data
    DETRAY_HOST
    explicit multi_axis(vecmem::memory_resource &resource) : m_data(resource) {}

    /// Constructor from data containers - move
    DETRAY_HOST_DEVICE
    multi_axis(vector_type<dindex_range> &&axes_data,
               vector_type<scalar_type> &&edges)
        : m_data(std::move(axes_data), std::move(edges)) {}

    /// Constructor from externally owned data containers.
    DETRAY_HOST_DEVICE
    multi_axis(const vector_type<dindex_range> *axes_data_ptr,
               const vector_type<scalar_type> *edges_ptr, dindex offset = 0)
        : m_data(axes_data_ptr, edges_ptr, offset) {}

    /// Device-side construction from a vecmem based view type
    template <typename axes_view_t,
              typename std::enable_if_t<
                  detray::detail::is_device_view_v<axes_view_t>, bool> = true>
    DETRAY_HOST_DEVICE multi_axis(const axes_view_t &view)
        : m_data(detray::detail::get<0>(view.m_view),
                 detray::detail::get<1>(view.m_view)) {}

    /// @returns the underlying axes storage. Either the container
    /// or a container pointer to a global collection - const
    DETRAY_HOST_DEVICE
    auto data() const -> const storage_type & { return m_data; }

    /// @returns the underlying axes storage. Either the container
    /// or a container pointer to a global collection - non-const for vecmem
    // TODO: Don't do
    DETRAY_HOST_DEVICE
    auto data() -> storage_type & { return m_data; }

    /// Build an axis object in place.
    ///
    /// @tparam index the position of the axis in the parameter pack. Also
    ///               determines which axes data are used to build the instance.
    ///
    /// @returns an axis object, corresponding to the index.
    template <std::size_t index>
    DETRAY_HOST_DEVICE typename label_matcher<axis_reg::to_id(index)>::type
    get_axis() const {
        return {m_data.axis_data(index), m_data.edges()};
    }

    /// Build an axis object in place.
    ///
    /// @tparam L label of the axis.
    ///
    /// @returns an axis object, corresponding to the label.
    template <n_axis::label L>
    DETRAY_HOST_DEVICE typename label_matcher<L>::type get_axis() const {
        return get_axis<axis_reg::to_index(L)>();
    }

    /// Build an axis object in place.
    ///
    /// @tparam axis_t type of the axis.
    ///
    /// @returns an axis object of the given type.
    template <typename axis_t>
    DETRAY_HOST_DEVICE axis_t get_axis() const {
        return get_axis<axis_t::shape_type::label>();
    }

    /// @returns the number of bins per axis
    DETRAY_HOST_DEVICE inline constexpr auto nbins() const -> multi_bin<Dim> {
        // Empty bin indices to be filled
        multi_bin<Dim> n_bins{{0, 0, 0}};
        // Get the number of bins for every axis
        (single_axis(get_axis<axis_ts>(), n_bins), ...);

        return n_bins;
    }

    /// Query the bin index for every coordinate of the given point on the axes.
    ///
    /// @tparam point_t the point in the local coordinate system that is spanned
    ///                 by the axes.
    ///
    /// @returns a multi bin that contains the resulting bin indices for
    ///          every axis in the corresponding entry (e.g. bin_x in entry 0)
    template <typename point_t>
    DETRAY_HOST_DEVICE multi_bin<Dim> bins(const point_t &p) const {
        // Empty bin indices to be filled
        multi_bin<Dim> bin_indices{};
        // Run the bin resolution for every axis in this multi-axis type
        (single_axis(get_axis<axis_ts>(), p, bin_indices), ...);

        return bin_indices;
    }

    /// @brief Get a bin range on every axis corresponding to the given point
    /// and the neighborhood around it.
    ///
    /// The neighborhood around the point can be defined in two ways:
    /// - scalar: neighborhood around the axis value,
    /// - index: neighborhood in #bins
    /// The resulting bin index range will contain all bins that belong to a
    /// given neighborhood around the lookup point.
    ///
    /// @tparam point_t the point in the local coordinate system that is spanned
    ///                 by the axes.
    /// @tparam neighbor_t the type of neighborhood defined on the axis around
    ///                    the point
    ///
    /// @returns a multi bin range that contains the resulting bin ranges for
    ///          every axis in the corresponding entry (e.g. rng_x in entry 0)
    template <typename point_t, typename neighbor_t>
    DETRAY_HOST_DEVICE multi_bin_range<Dim> bin_ranges(
        const point_t &p, const std::array<neighbor_t, 2> &nhood) const {
        // Empty bin ranges to be filled
        multi_bin_range<Dim> bin_ranges{};
        // Run the range resolution for every axis in this multi-axis type
        (single_axis(get_axis<axis_ts>(), p, nhood, bin_ranges), ...);

        return bin_ranges;
    }

    /// @returns a vecmem view on the axes data. Only allowed if it owning data.
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST auto get_data() -> view_type {
        return {detray::get_data(m_data.m_axes_data),
                detray::get_data(m_data.m_edges)};
    }

    /// @returns a vecmem const view on the axes data. Only allowed if it
    /// owning data.
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST auto get_data() const -> const_view_type {
        return {detray::get_data(m_data.m_axes_data),
                detray::get_data(m_data.m_edges)};
    }

    private:
    /// Get the number of bins for a single axis.
    ///
    /// @tparam axis_t defines the axis for the lookup (axis types are unique)
    ///
    /// @param ax the axis that performs the lookup
    /// @param n_bins the resulting bin numbers
    template <typename axis_t>
    DETRAY_HOST_DEVICE void single_axis(const axis_t &ax,
                                        multi_bin<Dim> &n_bins) const {
        // Get the index corresponding to the axis label (e.g. bin_x <=> 0)
        constexpr dindex loc_idx =
            axis_reg::to_index(axis_t::shape_type::label);
        n_bins[loc_idx] = ax.nbins();
    }

    /// Perform the bin lookup on a particular axis
    ///
    /// @tparam axis_t defines the axis for the lookup (axis types are unique)
    /// @tparam point_t is the point type in the axes local coordinates.
    ///
    /// @param ax the axis that performs the lookup
    /// @param p the point to be looked up on the axis
    /// @param bin_indices the multi-bin object that is filled with the results
    ///                    (axis index corresponds to entry of the multi-bin
    ///                    (e.g. binx <=> 0)
    template <typename axis_t, typename point_t>
    DETRAY_HOST_DEVICE void single_axis(const axis_t &ax, const point_t &p,
                                        multi_bin<Dim> &bin_indices) const {
        // Get the index corresponding to the axis label (e.g. bin_x <=> 0)
        constexpr dindex loc_idx =
            axis_reg::to_index(axis_t::shape_type::label);
        bin_indices.indices[loc_idx] = ax.bin(p[loc_idx]);
    }

    /// Perform the bin lookup on a particular axis within a given bin
    /// neighborhood
    ///
    /// @tparam axis_t defines the axis for the lookup (axis types are unique)
    /// @tparam point_t is the point type in the axes local coordinates.
    /// @tparam neighbor_t the type of neighborhood defined on the axis around
    ///                    the point
    ///
    /// @param ax the axis that performs the lookup
    /// @param p the point to be looked up on the axis
    /// @param nhood the neighborhood around the point for the range lookup
    /// @param bin_indices the multi-bin object that is filled with the results
    ///                    (axis index corresponds to entry of the multi-range
    ///                    (e.g. bin_rangex <=> 0)
    template <typename axis_t, typename point_t, typename neighbor_t>
    DETRAY_HOST_DEVICE void single_axis(
        const axis_t &ax, const point_t &p,
        const std::array<neighbor_t, 2> &nhood,
        multi_bin_range<Dim> &bin_ranges) const {
        // Get the index corresponding to the axis label (e.g. bin_range_x = 0)
        constexpr dindex loc_idx =
            axis_reg::to_index(axis_t::shape_type::label);
        bin_ranges.indices[loc_idx] = ax.range(p[loc_idx], nhood);
    }

    /// Internal storage of axes boundary ranges and bin edges
    storage_type m_data{};
};

}  // namespace n_axis

}  // namespace detray