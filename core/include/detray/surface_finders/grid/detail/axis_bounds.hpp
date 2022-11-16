/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/grid_axis.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// System include(s).
#include <cstddef>
#include <type_traits>

namespace detray::n_axis {

/// @brief Helper to tie two bin indices to a range.
///
/// @note Cannot use dindex_range for signed integer bin indices.
struct bin_range {
    /// lower bin in the range
    int lower{0};
    /// upper bin in the range
    int upper{0};

    /// Default constructor.
    constexpr bin_range() = default;

    /// Construction from two bin indices @param l (lower) @param u (upper).
    DETRAY_HOST_DEVICE
    constexpr bin_range(const int l, const int u) : lower{l}, upper{u} {}

    /// Implicit conversion from an @tparam array_t container of bin indices.
    template <template <typename, std::size_t> class array_t>
    DETRAY_HOST_DEVICE constexpr bin_range(const array_t<int, 2>& bin_array)
        : lower{bin_array[0]}, upper{bin_array[1]} {}

    /// Equality operator @returns true if both bin indices match.
    DETRAY_HOST_DEVICE
    constexpr bool operator==(const bin_range& rhs) const noexcept {
        return (lower == rhs.lower && upper == rhs.upper);
    }

    /// @returns the distance between the two bin indices.
    DETRAY_HOST_DEVICE
    constexpr int nbins() const { return upper - lower; }
};

/// @brief Describes the behaviour of an open axis.
///
/// The axis will be open, i.e. each underflow bin is mapped to 0 and each
/// overflow bin is mapped to #bins + 1: [0, #bins + 1]. Where 0 and #bins + 1
/// are the overflow bins.
///
/// @tparam axis_label the label of the axis, i.e. x, y, z or r.
template <n_axis::label axis_label>
struct open {

    static constexpr n_axis::label label = axis_label;
    static constexpr n_axis::bounds type = bounds::e_open;

    /// Access function to a single bin from a value v.
    ///
    /// @param ibin bin index to be mapped to axis bounds
    /// @param nbins is the total number of bins
    ///
    /// @returns an open axis bin index
    DETRAY_HOST_DEVICE
    constexpr auto map(const int ibin, const std::size_t nbins) const noexcept
        -> dindex {
        if (ibin <= 0) {
            // underflow bin
            return 0;
        } else if (ibin >= static_cast<int>(nbins)) {
            // overflow bin
            return static_cast<dindex>(nbins + 1);
        } else {
            // Shift the regular bins into the range [1, #bins]
            return static_cast<dindex>(ibin + 1);
        }
    }

    /// Access function to a range of bins.
    ///
    /// @param lbin the lower bin of the range.
    /// @param ubin is the upper bin of the range.
    /// @param nbins is the total number of bins
    ///
    /// @returns open bin range
    DETRAY_HOST_DEVICE
    constexpr auto map(const int lbin, const int ubin,
                       const std::size_t nbins) const noexcept -> dindex_range {

        dindex min_bin = (lbin >= 0) ? static_cast<dindex>(lbin + 1) : 0;
        dindex max_bin = (ubin < static_cast<int>(nbins))
                             ? static_cast<dindex>(ubin + 1)
                             : static_cast<dindex>(nbins + 1);

        return {min_bin, max_bin};
    }

    /// Access function to a range of bins - convenience function
    ///
    /// @param range signed range to be mapped to axis bounds
    /// @param nbins is the total number of bins
    ///
    /// @returns open bin range
    DETRAY_HOST_DEVICE
    constexpr auto map(const bin_range range,
                       const std::size_t nbins) const noexcept -> dindex_range {
        return map(range.lower, range.upper, nbins);
    }
};

/// @brief Describes the behaviour of a closed axis.
///
/// The axis will be closed, i.e. each underflow bin is mapped to 0 and each
/// overflow bin is mapped to #bins - 1. [0, #bins - 1]. Meaning, there are no
/// actual over- or underflow bins (they would be -1 and #bins).
///
/// @tparam axis_label the label of the axis, i.e. x, y, z or r.
template <n_axis::label axis_label>
struct closed {

    static constexpr n_axis::label label = axis_label;
    static constexpr bounds type = bounds::e_closed;

    /// Access function to a single bin from a value v
    ///
    /// @param ibin bin index to be mapped to axis bounds
    /// @param nbins is the total number of bins
    ///
    /// @returns a closed axis bin index
    DETRAY_HOST_DEVICE
    constexpr auto map(const int ibin, const std::size_t nbins) const noexcept
        -> dindex {
        if (ibin <= 0) {
            // underflow gets mapped onto axis bin 0
            return 0;
        } else if (ibin >= static_cast<int>(nbins)) {
            // overflow gets mapped onto axis bin #bins - 1
            return static_cast<dindex>(static_cast<int>(nbins) - 1);
        } else {
            return static_cast<dindex>(ibin);
        }
    }

    /// Access function to a range of bins
    ///
    /// @param lbin the lower bin of the range
    /// @param ubin is the upper bin of the range
    /// @param nbins is the total number of bins
    ///
    /// @returns closed bin range
    DETRAY_HOST_DEVICE
    inline dindex_range map(const int lbin, const int ubin,
                            const std::size_t nbins) const {

        dindex min_bin = (lbin > 0) ? static_cast<dindex>(lbin) : 0;
        dindex max_bin = (ubin >= static_cast<int>(nbins))
                             ? static_cast<dindex>(static_cast<int>(nbins) - 1)
                             : static_cast<dindex>(ubin);

        return {min_bin, max_bin};
    }

    /// Access function to a range of bins - convenience function
    ///
    /// @param range signed range to be mapped to axis bounds
    /// @param nbins is the total number of bins
    ///
    /// @returns closed bin range
    DETRAY_HOST_DEVICE
    inline dindex_range map(const bin_range range,
                            const std::size_t nbins) const {
        return map(range.lower, range.upper, nbins);
    }
};

/// @brief Describes the behaviour of a circular axis.
///
/// The axis will be periodic, i.e. underflow bins map into #bins - 1 and
/// overflow bins map into 0: so that [0, #bins - 1], with -1 = #bins - 1 and
/// #bins = 0.
template <n_axis::label axis_label = n_axis::label::e_phi>
struct circular {

    static constexpr n_axis::label label = axis_label;
    static constexpr bounds type = bounds::e_circular;

    /// Access function to a single bin from a value v
    ///
    /// @param ibin bin index to be mapped to axis bounds
    /// @param nbins is the total number of bins
    ///
    /// @returns a circular axis bin index
    DETRAY_HOST_DEVICE
    auto constexpr map(const int ibin, const std::size_t nbins) const noexcept
        -> dindex {
        if (ibin >= 0 and ibin < static_cast<int>(nbins)) {
            return static_cast<dindex>(ibin);
        } else {
            return static_cast<dindex>(wrap(ibin, nbins));
        }
    }

    /// Access function to a range of bins
    ///
    /// @param lbin the lower bin of the range
    /// @param ubin is the upper bin of the range
    /// @param nbins is the total number of bins
    ///
    /// The axis is circular: it @returns an ordered dindex_range: If the
    /// second range index is larger than the first, there has been a wraparound
    DETRAY_HOST_DEVICE
    auto constexpr map(const int lbin, const int ubin,
                       const std::size_t nbins) const noexcept -> dindex_range {
        dindex min_bin = static_cast<dindex>(wrap(lbin, nbins));
        dindex max_bin = static_cast<dindex>(wrap(ubin, nbins));
        return {min_bin, max_bin};
    }

    /// Access function to a range of bins - convenience function
    ///
    /// @param range signed range to be mapped to axis bounds
    /// @param nbins is the total number of bins
    ///
    /// The axis is circular: it @returns an ordered dindex_range: If the
    /// second range index is larger than the first, there has been a wraparound
    DETRAY_HOST_DEVICE
    auto constexpr map(const bin_range range,
                       const std::size_t nbins) const noexcept -> dindex_range {
        return map(range.lower, range.upper, nbins);
    }

    /// Wraps the bin index around for the periodic boundary condition
    ///
    /// @param ibin bin index to be mapped to axis bounds
    /// @param nbins is the total number of bins
    ///
    /// @return an index of a remapped bin
    DETRAY_HOST_DEVICE
    auto constexpr wrap(const int ibin, const std::size_t nbins) const -> int {
        const int bins = static_cast<int>(nbins);
        return (bins + (ibin % bins)) % bins;
    }
};

}  // namespace detray::n_axis