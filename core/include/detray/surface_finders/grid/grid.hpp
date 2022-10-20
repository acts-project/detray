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
#include "detray/surface_finders/grid/detail/grid_helpers.hpp"
#include "detray/surface_finders/grid/populator.hpp"
#include "detray/surface_finders/grid/serializer.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <cstddef>
#include <type_traits>

namespace detray {

/// @brief An N-dimensional grid for object storage.
///
/// @tparam multi_axis_t the types of grid axes
/// @tparam value_t type of data to be stored in the grid's bins.
/// @tparam serialzier_t how to serialize axis-local bin indices into global bin
///                      indices in the grid backend storage and vice versa.
/// @tparam populator_impl_t is a prescription what to do when a bin gets
///                          populated or read. It thus also defines the backend
///                          storage type/value type.
template <typename multi_axis_t, typename value_t,
          template <std::size_t> class serializer_t, typename populator_impl_t>
class grid {

    public:
    // Single value in a bin entry
    using value_type = value_t;
    template <std::size_t DIM>
    using serializer_type = serializer_t<DIM>;
    // Interface to the populator (detramines the bin entry and value type).
    using populator_type = populator<populator_impl_t>;
    using bin_type = typename populator_type::template bin_type<value_type>;
    /// The type of the multi-axis is tied to the type of the grid: a non-
    /// owning grid holds a non-owning multi-axis member.
    using axes_type = multi_axis_t;
    using container_types = typename axes_type::container_types;
    using local_frame = typename axes_type::local_frame_type;
    using scalar_type = typename axes_type::scalar_type;

    /// Grid dimension
    static constexpr std::size_t Dim = axes_type::Dim;
    static constexpr bool is_owning = axes_type::is_owning;

    /// How to define a neighborhood for this grid
    template <typename neighbor_t>
    using neighborhood_type =
        typename container_types::template array_type<neighbor_t, Dim>;

    /// Backend storage type for the grid
    using bin_storage_type =
        typename container_types::template vector_type<bin_type>;
    /// Vecmem based grid view type
    using view_type =
        dmulti_view<dvector_view<bin_type>, typename axes_type::view_type>;
    /// Vecmem based grid view type
    using const_view_type = dmulti_view<dvector_view<const bin_type>,
                                        typename axes_type::const_view_type>;
    /// Grid backend can be owning (single grid) or non-owning (grid collection)
    using storage_type =
        std::conditional_t<is_owning, detail::grid_data<bin_storage_type>,
                           detail::grid_view<bin_storage_type>>;
    template <bool owning>
    using type = grid<typename multi_axis_t::template type<owning>, value_t,
                      serializer_t, populator_impl_t>;

    /// Make grid default constructible: Empty grid with empty axis
    grid() = default;

    /// Create empty grid with empty axes from specific vecmem memory resource
    DETRAY_HOST
    explicit grid(vecmem::memory_resource &resource)
        : m_data(resource), m_axes(resource) {}

    /// Create grid with well defined @param axes and @param bins_data - move
    DETRAY_HOST_DEVICE
    grid(bin_storage_type &&bin_data, axes_type &&axes)
        : m_data(std::move(bin_data)), m_axes(std::move(axes)) {}

    /// Create grid from container pointers - non-owning (both grid and axes)
    DETRAY_HOST_DEVICE
    grid(bin_storage_type *bin_data_ptr, axes_type &axes,
         const dindex offset = 0)
        : m_data(bin_data_ptr, offset), m_axes(axes) {}

    /// Create grid from container pointers - non-owning (both grid and axes)
    DETRAY_HOST_DEVICE
    grid(bin_storage_type *bin_data_ptr, axes_type &&axes,
         const dindex offset = 0)
        : m_data(bin_data_ptr, offset), m_axes(axes) {}

    /// Device-side construction from a vecmem based view type
    template <typename grid_view_t,
              typename std::enable_if_t<detail::is_device_view_v<grid_view_t>,
                                        bool> = true>
    DETRAY_HOST_DEVICE grid(grid_view_t &view)
        : m_data(detray::detail::get<0>(view.m_view)),
          m_axes(detray::detail::get<1>(view.m_view)) {}

    /// @returns the underlying bin content storage. Either the container
    /// or a container pointer to a global collection - const
    DETRAY_HOST_DEVICE
    auto data() const -> const storage_type & { return m_data; }

    /// @returns the multi-axis used by the grid - const
    DETRAY_HOST_DEVICE
    auto axes() const -> const axes_type & { return m_axes; }

    /// @returns an axis object, corresponding to the label.
    template <std::size_t index>
    DETRAY_HOST_DEVICE inline constexpr auto get_axis() const {
        return m_axes.template get_axis<index>();
    }

    /// @returns an axis object, corresponding to the label.
    template <n_axis::label L>
    DETRAY_HOST_DEVICE inline constexpr auto get_axis() const {
        return m_axes.template get_axis<L>();
    }

    /// @returns an axis object of the given type.
    template <typename axis_t>
    DETRAY_HOST_DEVICE inline constexpr axis_t get_axis() const {
        return m_axes.template get_axis<axis_t>();
    }

    /// @returns the total number of bins in the grid
    DETRAY_HOST_DEVICE inline constexpr auto nbins() const -> std::size_t {
        const auto n_bins_per_axis = m_axes.nbins();
        std::size_t n_bins{1};
        for (std::size_t i{0}; i < Dim; ++i) {
            n_bins *= n_bins_per_axis[i];
        }
        return n_bins;
    }

    /// Transform a point in global cartesian coordinates to local coordinates
    ///
    /// @param trf the placement transform of the grid (e.g. from a volume or
    ///            a surface).
    /// @param p   the point in global coordinates
    ///
    /// @returns a point in the coordinate system that is spanned by the grid's
    /// axes.
    template <typename transform_t, typename point3_t, typename vector3_t>
    auto global_to_local(const transform_t &trf, const point3_t &p,
                         const vector3_t &d) const {
        return local_frame().global_to_local(trf, p, d);
    }

    /// @returns the iterable view of the bin content
    /// @{
    /// @param indices the single indices corresponding to a multi_bin
    template <typename... I, std::enable_if_t<sizeof...(I) == Dim, bool> = true>
    DETRAY_HOST_DEVICE auto at(I... indices) const {
        return at(n_axis::multi_bin<Dim>{{indices...}});
    }

    /// @param mbin the multi-index of bins over all axes
    DETRAY_HOST_DEVICE
    auto at(const n_axis::multi_bin<Dim> &mbin) const {
        return at(m_serializer(m_axes, mbin));
    }

    /// @param gbin the multi-index of bins over all axes
    DETRAY_HOST_DEVICE
    auto at(dindex gbin) const {
        return m_populator.view(*(data().bin_data()), gbin + data().offset());
    }
    /// @}

    /// Find the value of a single bin
    ///
    /// @param p is point in the local frame
    ///
    /// @return the iterable view of the bin content
    template <typename point_t,
              std::enable_if_t<std::is_class_v<point_t>, bool> = true>
    DETRAY_HOST_DEVICE auto search(const point_t &p) const {
        return at(m_axes.bins(p));
    }

    /// @brief Return a neighborhood of values from the grid
    ///
    /// The lookup is done with a neighborhood around the bin which contains the
    /// point
    ///
    /// @tparam sort is a directive whether return a sorted neighborhood
    ///
    /// @param p is point in the local frame
    /// @param nhood is the binned/scalar neighborhood
    ///
    /// @return the sequence of values
    template <typename point_t, typename neighbor_t, bool sort = false,
              std::enable_if_t<std::is_class_v<point_t>, bool> = true>
    DETRAY_HOST_DEVICE auto search(const point_t &p,
                                   neighborhood_type<neighbor_t> &nhood) const
        -> void {
        n_axis::multi_bin_range<Dim> bin_ranges = axes().bin_ranges(p, nhood);
        // Return iterable over bin values in the multi-range

        // Placeholder
    }

    /// Poupulate a bin with a single one of its corresponding values @param v
    /// @{
    /// @param mbin the multi bin index to be populated
    DETRAY_HOST_DEVICE auto populate(const n_axis::multi_bin<Dim> mbin,
                                     const value_type &v) -> void {
        populate(m_serializer(m_axes, mbin), v);
    }

    /// @param mbin the multi bin index to be populated
    DETRAY_HOST_DEVICE auto populate(const n_axis::multi_bin<Dim> mbin,
                                     value_type &&v) -> void {
        populate(m_serializer(m_axes, mbin), std::move(v));
    }

    /// @param gbin the global bin index to be populated
    DETRAY_HOST_DEVICE auto populate(const dindex gbin, const value_type &v)
        -> void {
        m_populator(*(m_data.bin_data()), gbin + m_data.offset(), v);
    }

    /// @param gbin the global bin index to be populated
    DETRAY_HOST_DEVICE auto populate(const dindex gbin, value_type &&v)
        -> void {
        m_populator(*(m_data.bin_data()), gbin + m_data.offset(), std::move(v));
    }

    /// @param p the point in local coordinates that defines the bin to be
    ///          populated
    template <typename point_t,
              std::enable_if_t<std::is_class_v<point_t>, bool> = true>
    DETRAY_HOST_DEVICE auto populate(const point_t &p, value_type &&v) -> void {
        populate(m_serializer(m_axes, m_axes.bins(p)), std::move(v));
    }

    /// @param p the point in local coordinates that defines the bin to be
    ///          populated
    template <typename point_t,
              std::enable_if_t<std::is_class_v<point_t>, bool> = true>
    DETRAY_HOST_DEVICE auto populate(const point_t &p, const value_type &v)
        -> void {
        populate(m_serializer(m_axes, m_axes.bins(p)), v);
    }
    /// @}

    static constexpr auto serializer() -> serializer_t<Dim> { return {}; }

    /// @returns view of a grid, including the grids mulit_axis. Also valid if
    /// the value type of the grid is cv qualified (then value_t propagates
    /// quialifiers) - non-const
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST auto get_data() -> view_type {
        return {detray::get_data(*(m_data.bin_data())),
                detray::get_data(m_axes)};
    }

    /// @returns view of a grid, including the grids mulit_axis. Also valid if
    /// the value type of the grid is cv qualified (then value_t propagates
    /// quialifiers) - const
    template <bool owner = is_owning, std::enable_if_t<owner, bool> = true>
    DETRAY_HOST auto get_data() const -> const_view_type {
        return {detray::get_data(*(data().bin_data())),
                detray::get_data(axes())};
    }

    private:
    /// Struct that contains the grid's data state
    storage_type m_data{};
    /// The axes of the grid
    axes_type m_axes{};
    /// How to write and fetch bin values from the backend storage
    populator_type m_populator{};
    /// Serialization/Deserialization of multi-bin indices to the global bin
    /// data storage
    serializer_t<Dim> m_serializer{};
};

}  // namespace detray
