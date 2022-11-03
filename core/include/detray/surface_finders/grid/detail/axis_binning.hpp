/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// System include(s).
#include <cstddef>

namespace detray::n_axis {

/// axis binning type names.
enum class binning {
    e_regular = 0,
    e_irregular = 1,
};

/// @brief A regular binning scheme.
///
/// The binning on the input parameter space is regular and therefore only needs
/// the span and the number of bins to match a value to a bin.
template <typename dcontainers = host_container_types,
          typename scalar_t = scalar>
struct regular {

    // Extract container types
    using scalar_type = scalar_t;
    using container_types = dcontainers;
    template <typename T>
    using vector_type = typename dcontainers::template vector_type<T>;
    template <typename T, std::size_t N>
    using array_type = typename dcontainers::template array_type<T, N>;

    static constexpr binning type = binning::e_regular;

    /// For the regular binning, the range contains the offset into the
    /// bin edges container and the total number of bins. This is because
    /// the number of entries in the bin edge container is always two and
    /// instead the number of bins has to be known.
    const dindex_range *m_edges_range{nullptr};
    /// Access to the bin edges
    const vector_type<scalar_t> *m_bin_edges{nullptr};

    /// Default constructor (no concrete memory acces)
    regular() = default;

    /// Constructor from an index range and bin edges - non-owning
    ///
    /// @param range range of bin boundary entries in an external storage
    /// @param edges lower edges for all bins in an external storage
    DETRAY_HOST_DEVICE
    regular(const dindex_range *range, const vector_type<scalar_t> *edges)
        : m_edges_range(range), m_bin_edges(edges) {}

    /// @returns the total number of bins, which for the regular axis is simply
    /// the second entry in the range
    DETRAY_HOST_DEVICE
    std::size_t nbins() const { return detail::get<1>(*m_edges_range); }

    /// Access function to a single bin from a value v
    ///
    /// @param v is the value for the bin search
    ///
    /// @note the floating point truncation to integer leaves
    /// int(-0.9) = int(0.9) = 0, which wrongfully maps overflow bins onto the
    /// axis range itself.
    /// Therefore, the values are shifted towards positive bin numbers first,
    /// and then shifted back, in order to emulate @c std::floor for this case.
    ///
    /// @returns the corresponding bin index
    DETRAY_HOST_DEVICE
    int bin(const scalar_t v) const {
        return static_cast<int>((v - span()[0]) / bin_width() + 1) - 1;
    }

    /// Access function to a range with binned neighborhood
    ///
    /// @param v is the value for the bin search
    /// @param nhood is the neighborhood bin index range (# neighboring bins)
    ///
    /// @returns the corresponding range of bin indices
    DETRAY_HOST_DEVICE
    array_type<int, 2> range(const scalar_t v,
                             const array_type<dindex, 2> &nhood) const {
        const int ibin{bin(v)};
        const int ibinmin{ibin - static_cast<int>(nhood[0])};
        const int ibinmax{ibin + static_cast<int>(nhood[1])};

        return {ibinmin, ibinmax};
    }

    /// Access function to a range with scalar neighborhood
    ///
    /// @param v is the value for the bin search
    /// @param nhood is the neighborhood value range (range on axis values)
    ///
    /// @returns the corresponding range of bin indices
    DETRAY_HOST_DEVICE
    array_type<int, 2> range(const scalar_t v,
                             const array_type<scalar_t, 2> &nhood) const {
        return {bin(v - nhood[0]), bin(v + nhood[1])};
    }

    /// @return the bin edges for a given @param ibin
    DETRAY_HOST_DEVICE
    array_type<scalar_t, 2> bin_edges(const dindex ibin) const {
        const scalar_t width{bin_width()};
        const scalar_t lower_edge{span()[0] + ibin * width};
        return {lower_edge, lower_edge + width};
    }

    /// @return the values of the edges of all bins - uses dynamic memory
    // TODO: return generator view instead to make it work properly in device
    DETRAY_HOST_DEVICE
    vector_type<scalar_t> bin_edges() const {
        // Output vector has to be constructed, because the edges are
        // calculated on the fly
        vector_type<scalar_t> edges;
        detail::call_reserve(edges, nbins());

        // Calculate bin edges from number of bins and axis span
        const array_type<scalar_t, 2> sp = span();
        const scalar_t step{bin_width()};

        for (dindex ib = 0; ib < nbins(); ++ib) {
            edges.push_back(sp[0] + ib * step);
        }

        return edges;
    }

    /// @return the bin width between any two bins.
    DETRAY_HOST_DEVICE
    scalar_t bin_width() const {
        // Get the binning information
        const dindex ibin{detail::get<0>(*m_edges_range)};
        const scalar_t min{(*m_bin_edges)[ibin]};
        const scalar_t max{(*m_bin_edges)[ibin + 1]};

        const scalar_t step_size{(max - min) / static_cast<scalar_t>(nbins())};

        return step_size;
    }

    /// @return the span of the binning (equivalent to the span of the axis:
    /// [min, max) )
    DETRAY_HOST_DEVICE
    array_type<scalar_t, 2> span() const {
        // Get the binning information
        const dindex ibin{detail::get<0>(*m_edges_range)};
        const scalar_t min{(*m_bin_edges)[ibin]};
        const scalar_t max{(*m_bin_edges)[ibin + 1]};

        return {min, max};
    }
};

/// @brief An irregular binning scheme.
///
/// The bin edges are irregular in the underlying parameter space.
/// Therefore, the correct bin index has to be searched for.
///
/// @note The bin search makes this type comparatively expensive. Only use when
/// absolutely needed.
template <typename dcontainers = host_container_types,
          typename scalar_t = scalar>
struct irregular {

    // Extract container types
    using scalar_type = scalar_t;
    using container_types = dcontainers;
    template <typename T>
    using vector_type = typename dcontainers::template vector_type<T>;
    template <typename T, std::size_t N>
    using array_type = typename dcontainers::template array_type<T, N>;

    static constexpr binning type = binning::e_irregular;

    /// Index range into the bin edges container
    const dindex_range *m_edges_range{nullptr};
    /// Access to the bin edges
    const vector_type<scalar_t> *m_bin_edges{nullptr};

    /// Default constructor (no concrete memory access)
    irregular() = default;

    /// Constructor from an index range and bin edges - non-owning
    ///
    /// @param range range of bin boundary entries in an external storage
    /// @param edges lower edges for all bins in an external storage
    DETRAY_HOST_DEVICE
    irregular(const dindex_range *range, const vector_type<scalar_t> *edges)
        : m_edges_range(range), m_bin_edges(edges) {}

    /// @returns the total number of bins
    DETRAY_HOST_DEVICE
    std::size_t nbins() const {
        return detail::get<1>(*m_edges_range) - detail::get<0>(*m_edges_range);
    }

    /// Access function to a single bin from a value v
    ///
    /// @param v is the value for the bin search
    ///
    /// @returns the corresponding bin index
    DETRAY_HOST_DEVICE
    int bin(const scalar_t v) const {
        auto bins_begin = m_bin_edges->begin() + detail::get<0>(*m_edges_range);
        auto bins_end = m_bin_edges->begin() + detail::get<1>(*m_edges_range);

        return static_cast<int>(detail::lower_bound(bins_begin, bins_end, v) -
                                bins_begin) -
               1;
    }

    /// Access function to a range with binned neighborhood
    ///
    /// @param v is the value for the bin search
    /// @param nhood is the neighborhood range (# neighboring bins)
    ///
    /// @returns the corresponding range of bin indices
    DETRAY_HOST_DEVICE
    array_type<int, 2> range(const scalar_t v,
                             const array_type<dindex, 2> &nhood) const {
        const int ibin{bin(v)};
        const int ibinmin{ibin - static_cast<int>(nhood[0])};
        const int ibinmax{ibin + static_cast<int>(nhood[1])};

        return {ibinmin, ibinmax};
    }

    /// Access function to a range with scalar neighborhood
    ///
    /// @param v is the value for the bin search
    /// @param nhood is the neighborhood range (range on axis values)
    ///
    /// @returns the corresponding range of bin indices
    DETRAY_HOST_DEVICE
    array_type<int, 2> range(const scalar_t v,
                             const array_type<scalar_t, 2> &nhood) const {
        const int nbin{bin(v - nhood[0])};
        const int pbin{bin(v + nhood[1])};
        return {nbin, pbin};
    }

    /// @return the bin edges for a given @param ibin
    DETRAY_HOST_DEVICE
    array_type<scalar_t, 2> bin_edges(const dindex ibin) const {
        // Offset into global container for this binning
        const dindex offset{detail::get<0>(*m_edges_range)};

        return {(*m_bin_edges)[offset + ibin],
                (*m_bin_edges)[offset + ibin + 1]};
    }

    /// @return the values of the edges of all bins - uses dynamic memory
    // TODO: return range view instead to make it work properly in device
    DETRAY_HOST_DEVICE
    vector_type<scalar_t> bin_edges() const {
        // Transcribe the subvector for this binning from the global storage
        vector_type<scalar_t> edges;
        detail::call_reserve(edges, nbins());

        edges.insert(edges->end(),
                     m_bin_edges->begin() + detail::get<0>(*m_edges_range),
                     m_bin_edges->begin() + detail::get<1>(*m_edges_range));

        return edges;
    }

    /// @return the bin width of a bin with index @param ibin.
    DETRAY_HOST_DEVICE
    scalar_t bin_width(const dindex ibin) const {

        // Get the binning information
        const array_type<scalar_t, 2> edges = bin_edges(ibin);

        return edges[1] - edges[0];
    }

    /// @return the span of the binning (equivalent to the span of the axis:
    /// [min, max) )
    DETRAY_HOST_DEVICE
    array_type<scalar_t, 2> span() const {
        // Get the binning information
        const scalar_t min{(*m_bin_edges)[detail::get<0>(*m_edges_range)]};
        const scalar_t max{(*m_bin_edges)[detail::get<1>(*m_edges_range)]};

        return {min, max};
    }
};

}  // namespace detray::n_axis
