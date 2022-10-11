/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// System include(s).
#include <cstddef>
#include <type_traits>

namespace detray {

namespace n_axis {

/// axis shape names.
enum class shape {
    e_open = 0,
    e_closed = 1,
    e_circular = 1,
};

/// axis coordinate names. Used to get a specific axes from an axes collection.
enum class label {
    e_x = 0,
    e_y = 1,
    e_z = 2,
    e_r = 0,
    e_phi = 1,
    e_rphi = 0,
    e_cyl_z = 1,
};

/// @brief Helper to tie two bin indices to a range.
///
/// @note Cannot use dindex_range for signed integer bin indices.
struct bin_range {
    /// lower bin in the range
    int lower{0};
    /// upper bin in the range
    int upper{0};

    /// Default constructor.
    bin_range() = default;

    /// Construction from two bin indices @param l (lower) @param u (upper).
    DETRAY_HOST_DEVICE
    bin_range(const int l, const int u) : lower{l}, upper{u} {}

    /// Implicit conversion from an @tparam array_t container of bin indices.
    template <template <typename, std::size_t> class array_t>
    DETRAY_HOST_DEVICE bin_range(const array_t<int, 2>& bin_array)
        : lower{bin_array[0]}, upper{bin_array[1]} {}

    /// Equality operator @returns true if both bin indices match.
    DETRAY_HOST_DEVICE
    inline bool operator==(const bin_range& rhs) const {
        return (lower == rhs.lower && upper == rhs.upper);
    }

    /// @returns the distance between the two bin indices.
    DETRAY_HOST_DEVICE
    int nbins() const { return upper - lower; }
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
    static constexpr n_axis::shape type = shape::e_open;

    /// Access function to a single bin from a value v.
    ///
    /// @param ibin bin index to be mapped to axis shape
    /// @param nbins is the total number of bins
    ///
    /// @returns an open axis bin index
    DETRAY_HOST_DEVICE inline dindex map(const int ibin,
                                         const std::size_t nbins) const {
        if (ibin <= 0) {
            return 0;
        } else if (ibin > static_cast<int>(nbins)) {
            return static_cast<dindex>(nbins + 1);
        } else {
            return static_cast<dindex>(ibin);
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
    inline dindex_range map(const int lbin, const int ubin,
                            const std::size_t nbins) const {

        dindex min_bin = (lbin > 0) ? static_cast<dindex>(lbin) : 0;
        dindex max_bin = (ubin <= static_cast<int>(nbins))
                             ? static_cast<dindex>(ubin)
                             : static_cast<dindex>(nbins + 1);

        return {min_bin, max_bin};
    }

    /// Access function to a range of bins - convenience function
    ///
    /// @param range signed range to be mapped to axis shape
    /// @param nbins is the total number of bins
    ///
    /// @returns open bin range
    DETRAY_HOST_DEVICE
    inline dindex_range map(const bin_range range,
                            const std::size_t nbins) const {
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
    static constexpr shape type = shape::e_closed;

    /// Access function to a single bin from a value v
    ///
    /// @param ibin bin index to be mapped to axis shape
    /// @param nbins is the total number of bins
    ///
    /// @returns a closed axis bin index
    DETRAY_HOST_DEVICE inline dindex map(const int ibin,
                                         const std::size_t nbins) const {
        if (ibin <= 0) {
            return 0;
        } else if (ibin > static_cast<int>(nbins) - 1) {
            return static_cast<dindex>(static_cast<int>(nbins) - 1);
        } else {
            return static_cast<dindex>(ibin - 1);
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

        dindex min_bin = (lbin > 0) ? static_cast<dindex>(lbin - 1) : 0;
        dindex max_bin = (ubin > static_cast<int>(nbins) - 1)
                             ? static_cast<dindex>(static_cast<int>(nbins) - 1)
                             : static_cast<dindex>(ubin - 1);

        return {min_bin, max_bin};
    }

    /// Access function to a range of bins - convenience function
    ///
    /// @param range signed range to be mapped to axis shape
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
/// The axis will be periodic, i.e. underflow and underflow bins map into 0:
/// 0 = #bins, so that [0, #bins - 1].
// TODO: rename to circular, once there is no name clash
template <n_axis::label axis_label = n_axis::label::e_phi>
struct circular {

    static constexpr n_axis::label label = axis_label;
    static constexpr shape type = shape::e_circular;

    /// Access function to a single bin from a value v
    ///
    /// @param ibin bin index to be mapped to axis shape
    /// @param nbins is the total number of bins
    ///
    /// @returns a circular axis bin index
    DETRAY_HOST_DEVICE inline dindex map(const int ibin,
                                         const std::size_t nbins) const {
        if (ibin > 0 and ibin < static_cast<int>(nbins)) {
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
    inline dindex_range map(const int lbin, const int ubin,
                            const std::size_t nbins) const {
        dindex min_bin = static_cast<dindex>(wrap(lbin, nbins));
        dindex max_bin = static_cast<dindex>(wrap(ubin, nbins));
        return {min_bin, max_bin};
    }

    /// Access function to a range of bins - convenience function
    ///
    /// @param range signed range to be mapped to axis shape
    /// @param nbins is the total number of bins
    ///
    /// The axis is circular: it @returns an ordered dindex_range: If the
    /// second range index is larger than the first, there has been a wraparound
    DETRAY_HOST_DEVICE
    inline dindex_range map(const bin_range range,
                            const std::size_t nbins) const {
        return map(range.lower, range.upper, nbins);
    }

    /// Wraps the bin index around for the periodic boundary condition
    ///
    /// @param ibin bin index to be mapped to axis shape
    /// @param nbins is the total number of bins
    ///
    /// @return an index of a remapped bin
    DETRAY_HOST_DEVICE inline int wrap(const int ibin,
                                       const std::size_t nbins) const {
        const int bins = static_cast<int>(nbins);
        return (bins + ((ibin - 1) % bins)) % bins;
    }
};

/// Determine axis shape as either 'open' or 'closed' for non-circular axes
template <n_axis::shape s, n_axis::label axis_label>
using shape_t =
    std::conditional_t<s == n_axis::shape::e_open, n_axis::open<axis_label>,
                       n_axis::closed<axis_label>>;

}  // namespace n_axis

}  // namespace detray