/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <algorithm>
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/memory_resource.hpp>

#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/indexing.hpp"

namespace detray {

namespace axis {
/** A regular closed axis.
 *
 * The axis is closed, i.e. each underflow bin is mapped to 0
 * and henceforth each overflow bin is mapped to bins-1
 */
template <template <typename, unsigned int> class array_type = darray,
          template <typename...> class vector_type = dvector>
struct regular {
    dindex n_bins;
    scalar min;
    scalar max;

    /* dummy boundary which is not used */
    vector_type<scalar> boundaries;

    /** Defualt constructor **/
    regular() = default;

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    regular(dindex axis_bins, scalar axis_min, scalar axis_max,
            vecmem::memory_resource &resource)
        : n_bins(axis_bins),
          min(axis_min),
          max(axis_max),
          boundaries(&resource) {}

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    regular(const regular &axis, vecmem::memory_resource &resource)
        : n_bins(axis.n_bins),
          min(axis.min),
          max(axis.max),
          boundaries(axis.boundaries, &resource) {}

    /** Constructor with axis_data **/
    template <
        typename axis_data_t,
        std::enable_if_t<!std::is_same_v<regular, axis_data_t>, bool> = true>
    DETRAY_HOST_DEVICE regular(axis_data_t &axis_data)
        : n_bins(axis_data.n_bins),
          min(axis_data.min),
          max(axis_data.max),
          boundaries(axis_data.boundaries) {}

    static constexpr unsigned int axis_identifier = 0;

    /** Return the number of bins */
    DETRAY_HOST_DEVICE
    dindex bins() const { return n_bins; }

    /** Access function to a single bin from a value v
     *
     * @param v is the value for the bin search
     *
     * As the axis is closed it @returns a dindex type
     **/
    DETRAY_HOST_DEVICE
    dindex bin(scalar v) const {
        int ibin = static_cast<int>((v - min) / (max - min) * n_bins);
        if (ibin >= 0 and ibin < static_cast<int>(n_bins)) {
            return static_cast<dindex>(ibin);
        } else {
            if (ibin < 0) {
                return 0;
            } else {
                return static_cast<dindex>(n_bins - 1);
            }
        }
    }

    /** Access function to a range with binned neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (in #bins)
     *
     * As the axis is closed it @returns a dindex_range
     **/
    DETRAY_HOST_DEVICE
    dindex_range range(scalar v,
                       const array_type<dindex, 2> &nhood = {0u, 0u}) const {

        int ibin = static_cast<int>((v - min) / (max - min) * n_bins);
        int ibinmin = ibin - static_cast<int>(nhood[0]);
        int ibinmax = ibin + static_cast<int>(nhood[1]);
        dindex min_bin = (ibinmin >= 0) ? static_cast<dindex>(ibinmin) : 0;
        dindex max_bin = (ibinmax < static_cast<int>(n_bins))
                             ? static_cast<dindex>(ibinmax)
                             : static_cast<dindex>(n_bins - 1);
        return {min_bin, max_bin};
    }

    /** Access function to a range with scalar neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (scalar)
     *
     * As the axis is closed it @returns a dindex_range
     **/
    DETRAY_HOST_DEVICE
    dindex_range range(scalar v, const array_type<scalar, 2> &nhood) const {
        int nbin =
            static_cast<int>((v - nhood[0] - min) / (max - min) * n_bins);
        int pbin =
            static_cast<int>((v + nhood[1] - min) / (max - min) * n_bins);
        dindex min_bin = (nbin >= 0) ? static_cast<dindex>(nbin) : 0;
        dindex max_bin = (pbin < static_cast<int>(n_bins))
                             ? static_cast<dindex>(pbin)
                             : static_cast<dindex>(n_bins - 1);
        return {min_bin, max_bin};
    }

    /** Access function to a zone with binned neighborhood
     *
     *
     * @tparam neighbor_t is the neighborhood size
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (in #bins)
     *
     * As the axis is closed it @returns a dindex_sequence
     **/
    template <typename neighbor_t>
    DETRAY_HOST_DEVICE dindex_sequence
    zone_t(scalar v, const array_type<neighbor_t, 2> &nhood) const {
        dindex_range nh_range = range(v, nhood);
        dindex_sequence sequence(static_cast<dindex_sequence::size_type>(
                                     nh_range[1] - nh_range[0] + 1),
                                 nh_range[0]);
        dindex m = 0;
        std::for_each(sequence.begin(), sequence.end(),
                      [&](auto &n) { n += m++; });
        return sequence;
    }

    /** Access function to a zone with binned neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (#bins)
     *
     * As the axis is closed it @returns a dindex_sequence
     **/
    DETRAY_HOST_DEVICE
    dindex_sequence zone(scalar v, const array_type<dindex, 2> &nhood) const {
        return zone_t<dindex>(v, nhood);
    }

    /** Access function to a zone with scalar neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (scalar)
     *
     * As the axis is closed it @returns a dindex_sequence
     **/
    DETRAY_HOST_DEVICE
    dindex_sequence zone(scalar v, const array_type<scalar, 2> &nhood) const {
        return zone_t<scalar>(v, nhood);
    }

    /** @return the bin boundaries for a given @param ibin */
    DETRAY_HOST_DEVICE
    array_type<scalar, 2> borders(dindex ibin) const {
        scalar step = (max - min) / n_bins;
        return {min + ibin * step, min + (ibin + 1) * step};
    }

    /** @return the values of the borders */
    DETRAY_HOST_DEVICE
    vector_type<scalar> all_borders() const {
        vector_type<scalar> borders;
        borders.reserve(n_bins + 1);
        scalar step = (max - min) / n_bins;
        for (dindex ib = 0; ib < n_bins + 1; ++ib) {
            borders.push_back(min + ib * step);
        }
        return borders;
    }

    /** @return the axis span [min, max) */
    DETRAY_HOST_DEVICE
    array_type<scalar, 2> span() const { return {min, max}; }
};

/** A regular circular axis.
 *
 * The axis is circular, i.e. the underflow bins map into the circular sequence
 */
template <template <typename, unsigned int> class array_type = darray,
          template <typename...> class vector_type = dvector>
struct circular {

    dindex n_bins;
    scalar min;
    scalar max;

    /* dummy boundary which is not used */
    vector_type<scalar> boundaries;

    static constexpr unsigned int axis_identifier = 1;

    /** Defualt constructor **/
    circular() = default;

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    circular(dindex axis_bins, scalar axis_min, scalar axis_max,
             vecmem::memory_resource &resource)
        : n_bins(axis_bins),
          min(axis_min),
          max(axis_max),
          boundaries(&resource) {}

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    circular(const circular &axis, vecmem::memory_resource &resource)
        : n_bins(axis.n_bins),
          min(axis.min),
          max(axis.max),
          boundaries(axis.boundaries, &resource) {}

    /** Constructor with axis_data **/
    template <
        typename axis_data_t,
        std::enable_if_t<!std::is_same_v<circular, axis_data_t>, bool> = true>
    DETRAY_HOST_DEVICE circular(axis_data_t &axis_data)
        : n_bins(axis_data.n_bins),
          min(axis_data.min),
          max(axis_data.max),
          boundaries(axis_data.boundaries) {}

    /** Return the number of bins */
    DETRAY_HOST_DEVICE
    dindex bins() const { return n_bins; }

    /** Access function to a single bin from a value v
     *
     * @param v is the value for the bin search
     *
     * As the axis is closed it @returns a dindex type
     **/
    DETRAY_HOST_DEVICE
    dindex bin(scalar v) const {
        int ibin = static_cast<int>((v - min) / (max - min) * n_bins);
        if (ibin >= 0 and ibin < static_cast<int>(n_bins)) {
            return static_cast<dindex>(ibin);
        } else {
            if (ibin < 0) {
                return static_cast<dindex>(n_bins + ibin);
            } else {
                return static_cast<dindex>(ibin - n_bins);
            }
        }
    }

    /** Access function to a range with binned neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the  neighborhood range (in #bins)
     *
     * As the axis is circular it @returns a dindex_range
     **/
    DETRAY_HOST_DEVICE
    dindex_range range(scalar v,
                       const array_type<dindex, 2> nhood = {0u, 0u}) const {
        dindex gbin = bin(v);
        dindex min_bin = remap(gbin, -static_cast<int>(nhood[0]));
        dindex max_bin = remap(gbin, static_cast<int>(nhood[1]));
        return {min_bin, max_bin};
    }

    /** Access function to a range with scalar neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (scalar)
     *
     * As the axis is circular it @returns a dindex_range
     **/
    DETRAY_HOST_DEVICE
    dindex_range range(scalar v, const array_type<scalar, 2> &nhood) const {
        dindex nbin = bin(v - nhood[0]);
        dindex pbin = bin(v + nhood[1]);
        return {nbin, pbin};
    }

    /** Access function to a zone with binned/scalar neighborhood
     *
     * @tparam neighbor_t is the neighborhood size
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (in #bins/scalar)
     *
     * As the axis is closed it @returns a dindex_sequence
     **/
    template <typename neighbor_t>
    DETRAY_HOST_DEVICE dindex_sequence
    zone_t(scalar v, const array_type<neighbor_t, 2> &nhood) const {
        dindex_range nh_range = range(v, nhood);
        if (nh_range[0] < nh_range[1]) {
            dindex_sequence sequence(static_cast<dindex_sequence::size_type>(
                                         nh_range[1] - nh_range[0] + 1),
                                     nh_range[0]);
            dindex m = 0;
            std::for_each(sequence.begin(), sequence.end(),
                          [&](auto &n) { n += m++; });
            return sequence;
        }
        dindex vl = static_cast<dindex>(n_bins - nh_range[0] + nh_range[1] + 1);
        dindex mi = 0;
        dindex mo = 0;
        dindex_sequence sequence(static_cast<dindex_sequence::size_type>(vl),
                                 nh_range[0]);
        std::for_each(sequence.begin(), sequence.end(), [&](auto &n) {
            n += mi++;
            if (n > n_bins - 1) {
                n = mo++;
            }
        });
        return sequence;
    }

    /** Access function to a zone with binned neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (#bins)
     *
     * As the axis is closed it @returns a dindex_sequence
     **/
    DETRAY_HOST_DEVICE
    dindex_sequence zone(scalar v, const array_type<dindex, 2> &nhood) const {
        return zone_t<dindex>(v, nhood);
    }

    /** Access function to a zone with scalar neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (scalar)
     *
     * As the axis is closed it @returns a dindex_sequence
     **/
    DETRAY_HOST_DEVICE
    dindex_sequence zone(scalar v, const array_type<scalar, 2> &nhood) const {
        return zone_t<scalar>(v, nhood);
    }

    /** Helper function to remap onto a circular range
     *
     * @param ibin is the optional binning value
     * @param shood is the sided neighbour hood
     *
     * @return an index, remapped bin
     **/
    DETRAY_HOST_DEVICE
    dindex remap(dindex ibin, int shood) const {
        int opt_bin = static_cast<int>(ibin) + shood;
        if (opt_bin >= 0 and opt_bin < static_cast<int>(n_bins)) {
            return static_cast<dindex>(opt_bin);
        }
        if (opt_bin < 0) {
            return static_cast<dindex>(n_bins + opt_bin);
        }
        return static_cast<dindex>(opt_bin - n_bins);
    }

    /** @return the bin boundaries for a given @param ibin */
    DETRAY_HOST_DEVICE
    array_type<scalar, 2> borders(dindex ibin) const {
        scalar step = (max - min) / n_bins;
        return {min + ibin * step, min + (ibin + 1) * step};
    }

    /** @return the values of the borders */
    DETRAY_HOST_DEVICE
    vector_type<scalar> all_borders() const {
        vector_type<scalar> borders;
        borders.reserve(n_bins + 1);
        scalar step = (max - min) / n_bins;
        for (dindex ib = 0; ib < n_bins + 1; ++ib) {
            borders.push_back(min + ib * step);
        }
        return borders;
    }

    /** @return the range  */
    DETRAY_HOST_DEVICE
    array_type<scalar, 2> span() const { return {min, max}; }
};

/** An iregular circular axis.
 *
 * The axis is closed, i.e. the underflow is mapped into the first,
 * the overflow is mapped into the last.
 */
template <template <typename, unsigned int> class array_type = darray,
          template <typename...> class vector_type = dvector>
struct irregular {

    /* dummy bin size, min and max */
    dindex n_bins;
    scalar min;
    scalar max;

    vector_type<scalar> boundaries;

    /** Defualt constructor **/
    irregular() = default;

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST irregular(vecmem::memory_resource &resource)
        : boundaries(&resource) {}

    /** Constructor with vecmem memory resource - rvalue **/
    DETRAY_HOST irregular(vector_type<scalar> &&bins,
                          vecmem::memory_resource &resource)
        : n_bins(bins.size() - 1),
          min(bins[0]),
          max(bins[n_bins]),
          boundaries(bins, &resource) {}

    /** Constructor with vecmem memory resource - lvalue **/
    DETRAY_HOST irregular(vector_type<scalar> &bins,
                          vecmem::memory_resource &resource)
        : n_bins(bins.size() - 1),
          min(bins[0]),
          max(bins[n_bins]),
          boundaries(bins, &resource) {}

    /** Constructor with vecmem memory resource **/
    DETRAY_HOST
    irregular(const irregular &axis, vecmem::memory_resource &resource)
        : n_bins(axis.n_bins),
          min(axis.min),
          max(axis.max),
          boundaries(axis.boundaries, &resource) {}

    /** Constructor with axis_data **/
    template <
        typename axis_data_t,
        std::enable_if_t<!std::is_same_v<irregular, axis_data_t>, bool> = true>
    DETRAY_HOST_DEVICE irregular(axis_data_t &axis_data)
        : n_bins(axis_data.n_bins),
          min(axis_data.min),
          max(axis_data.max),
          boundaries(axis_data.boundaries) {}

    static constexpr unsigned int axis_identifier = 2;

    /** Return the number of bins */
    DETRAY_HOST_DEVICE
    dindex bins() const { return static_cast<dindex>(boundaries.size() - 1); }

    /** Access function to a single bin from a value v
     *
     * @param v is the value for the bin search
     *
     * As the axis is closed it @returns a dindex type
     **/
    DETRAY_HOST_DEVICE
    dindex bin(scalar v) const {
        int ibin = static_cast<int>(
            std::lower_bound(boundaries.begin(), boundaries.end(), v) -
            boundaries.begin());
        if (ibin > 0 and ibin < static_cast<int>(boundaries.size())) {
            return static_cast<dindex>(--ibin);
        } else {
            if (ibin == 0) {
                return static_cast<dindex>(ibin);
            } else {
                return static_cast<dindex>(boundaries.size() - 2);
            }
        }
    }

    /** Access function to a range with binned neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (#bins)
     *
     * As the axis is closed it @returns a dindex_range
     **/
    DETRAY_HOST_DEVICE
    dindex_range range(scalar v,
                       const array_type<dindex, 2> &nhood = {0u, 0u}) const {

        dindex ibin = bin(v);
        int bins = boundaries.size() - 1;
        int ibinmin = ibin - static_cast<int>(nhood[0]);
        int ibinmax = ibin + static_cast<int>(nhood[1]);
        dindex min_bin = (ibinmin >= 0) ? static_cast<dindex>(ibinmin) : 0;
        dindex max_bin = (ibinmax < static_cast<int>(bins))
                             ? static_cast<dindex>(ibinmax)
                             : static_cast<dindex>(bins - 1);
        return {min_bin, max_bin};
    }

    /** Access function to a range with scalar neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (scalar)
     *
     * As the axis is closed it @returns a dindex_range
     **/
    DETRAY_HOST_DEVICE
    dindex_range range(scalar v, const array_type<scalar, 2> &nhood) const {
        dindex nbin = bin(v - nhood[0]);
        dindex pbin = bin(v + nhood[1]);
        return {nbin, pbin};
    }

    /** Access function to a zone with binned/scalar neighborhood
     *
     * @tparam neighbor_t is the neighborhood type
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (binned/scalar)
     *
     * As the axis is closed it @returns a dindex_sequence
     **/
    template <typename neighbor_t>
    DETRAY_HOST_DEVICE dindex_sequence
    zone_t(scalar v, const array_type<neighbor_t, 2> nhood) const {
        dindex_range nh_range = range(v, nhood);
        dindex_sequence sequence(static_cast<dindex_sequence::size_type>(
                                     nh_range[1] - nh_range[0] + 1),
                                 nh_range[0]);
        dindex m = 0;
        std::for_each(sequence.begin(), sequence.end(),
                      [&](auto &n) { n += m++; });
        return sequence;
    }

    /** Access function to a zone with binned neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (#bins)
     *
     * As the axis is closed it @returns a dindex_sequence
     **/
    DETRAY_HOST_DEVICE
    dindex_sequence zone(scalar v,
                         const array_type<dindex, 2> &nhood = {0, 0}) const {
        return zone_t<dindex>(v, nhood);
    }

    /** Access function to a zone with scalar neighborhood
     *
     * @param v is the value for the bin search
     * @param nhood is the neighborhood range (scalar)
     *
     * As the axis is closed it @returns a dindex_sequence
     **/
    DETRAY_HOST_DEVICE
    dindex_sequence zone(scalar v, const array_type<scalar, 2> &nhood) const {
        return zone_t<scalar>(v, nhood);
    }

    DETRAY_HOST_DEVICE
    /** @return the bin boundaries for a given @param ibin */
    array_type<scalar, 2> borders(dindex ibin) const {
        return {boundaries[ibin], boundaries[ibin + 1]};
    }

    /** @return the values of the borders of all bins */
    DETRAY_HOST
    vector_type<scalar> all_borders() const { return boundaries; }

    /** @return the range  */
    DETRAY_HOST_DEVICE
    array_type<scalar, 2> span() const {
        return {boundaries[0], boundaries[boundaries.size() - 1]};
    }
};

}  // namespace axis

/**
 * static implementation of axis data for device
 */
template <typename axis_type>
struct axis_data {

    axis_data(axis_type &axis)
        : n_bins(axis.n_bins),
          min(axis.min),
          max(axis.max),
          boundaries(vecmem::get_data(axis.boundaries)) {}

    dindex n_bins;
    scalar min;
    scalar max;
    vecmem::data::vector_view<scalar> boundaries;
};

/**
 * standalone function to get axis_data
 */
template <template <template <typename, unsigned int> class,
                    template <typename...> class>
          class axis_type,
          template <typename, unsigned int> class array_type,
          template <typename...> class vector_type>
inline axis_data<axis_type<array_type, vector_type>> get_data(
    axis_type<array_type, vector_type> &axis) {
    return axis;
}

}  // namespace detray
