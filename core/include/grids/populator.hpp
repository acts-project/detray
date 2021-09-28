/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <algorithm>
#include <limits>
#include <vecmem/containers/data/jagged_vector_buffer.hpp>
#include <vecmem/containers/data/jagged_vector_view.hpp>
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

#include "definitions/invalid_values.hpp"
#include "definitions/cuda_qualifiers.hpp"
#include "utils/indexing.hpp"
namespace detray {
/** A replace populator that swaps whatever current value in the
 * bin with the new one.
 *
 * @tparam value_type the type of a single stored object
 *
 * @note bare_value and store_value are identicial in this case
 **/
template <typename value_type = dindex,
          template <typename> class vector_type = dvector>
struct replace_populator {
    DETRAY_HOST_DEVICE
    replace_populator(const value_type invalid = invalid_value<value_type>())
        : kInvalid(invalid) {}

    value_type kInvalid;

    using bare_value = value_type;
    using store_value = value_type;
    using serialized_storage = vector_type<store_value>;

    using vector_view_t = vecmem::data::vector_view<store_value>;
    using vector_buffer_t = vecmem::data::vector_buffer<store_value>;
    using buffer_size_t = typename vector_view_t::size_type;

    /** Swap the stored value with a new bare value
     *
     * @param stored the stored value for the population
     * @param bvalue the new value to be added
     **/
    DETRAY_HOST_DEVICE
    void operator()(store_value &stored, bare_value &&bvalue) const {
        stored = std::move(bvalue);
    }

    /** Create a sequence of bare values, independent of the store_value.
     *
     * @param stored the stored value
     *
     * @return a sequence of bare values
     */
    DETRAY_HOST_DEVICE
    vector_type<bare_value> sequence(store_value &stored) const {
        if (stored != kInvalid) {
            return {stored};
        }
        return {};
    }

    /** Shift operation for unified memory block
     *
     * @param stored the stored value
     * @param offset is the shift offset
     *
     **/
    DETRAY_HOST_DEVICE
    void shift(store_value &stored, const bare_value &offset) const {
        stored += offset;
    }

    /** Return an initialized bin value
     */
    DETRAY_HOST_DEVICE
    store_value init() const { return kInvalid; }

    /** Return a vector view
     **/
    DETRAY_HOST
    static vector_view_t get_data(serialized_storage &data,
                                  vecmem::memory_resource &resource) {
        return vecmem::get_data(data);
    }
};

/** A complete populator that adds values to the internal
 * store array until it is completed, ignored afterwards.
 *
 * @tparam kDIM the dimension of the underlying stored array
 * @tparam kSORT a sorting flag
 * @tparam value_type the type of a single stored object
 * @tparam kInvalid the chosen invalid type
 *
 * @note bare_value and store_value are different in this case
 **/
template <unsigned int kDIM, bool kSORT = false, typename value_type = dindex,
          template <typename, unsigned int> class array_type = darray,
          template <typename> class vector_type = dvector>
struct complete_populator {
    DETRAY_HOST_DEVICE
    complete_populator(const value_type invalid = invalid_value<value_type>())
        : kInvalid(invalid) {}

    value_type kInvalid;

    using bare_value = value_type;
    using store_value = array_type<bare_value, kDIM>;
    using serialized_storage = vector_type<store_value>;

    using vector_view_t = vecmem::data::vector_view<store_value>;
    using vector_buffer_t = vecmem::data::vector_buffer<store_value>;
    using buffer_size_t = typename vector_view_t::size_type;

    /** Complete the stored value with a new bare value - for host
     *
     * @param stored the stored value for the population
     * @param bvalue the new value to be added
     **/
    DETRAY_HOST_DEVICE
    void operator()(store_value &stored, bare_value &&bvalue) const {
        for (auto &val : stored) {
            if (val == kInvalid) {
                val = bvalue;
                break;
            }
        }
        // no sort function in the cuda device
        // maybe can use thrust sort function
#if !defined(__CUDACC__)
        if (kSORT) {
            std::sort(stored.begin(), stored.end());
        }
#endif
    }

    /** Create a sequence of bare values, independent of the store_value.
     *
     * @param stored the stored value
     *
     * @return a sequence of bare values, @note it will ignore invalid entries
     */
    DETRAY_HOST_DEVICE
    vector_type<bare_value> sequence(store_value &stored) const {
        vector_type<bare_value> s;
        s.reserve(kDIM);
        for (const auto &val : stored) {
            if (val != kInvalid) {
                s.push_back(val);
            }
        }
        return s;
    }

    /** Shift operation for unified memory block
     *
     * @param stored the stored value
     * @param offset is the shift offset
     *
     **/
    DETRAY_HOST_DEVICE
    void shift(store_value &stored, const bare_value &offset) const {
        std::for_each(stored.begin(), stored.end(),
                      [&](auto &d) { d += offset; });
    }

    /** Return an initialized bin value
     **/
    DETRAY_HOST_DEVICE
    store_value init() const {

        store_value init_bin;
        for (auto &val : init_bin) {
            val = kInvalid;
        }
        return init_bin;
    }

    /** Return a vector view
     **/
    DETRAY_HOST
    static vector_view_t get_data(serialized_storage &data,
                                  vecmem::memory_resource &resource) {
        return vecmem::get_data(data);
    }
};

/** An attach populator that adds the new value to the
 *
 * @tparam kSORT the sorting directive
 * @tparam value_type the type of a single stored object
 *
 * @note bare_value and store_value are identicial in this case
 **/
template <bool kSORT = false, typename value_type = dindex,
          template <typename> class vector_type = dvector,
          template <typename> class jagged_vector_type = djagged_vector>
struct attach_populator {
    DETRAY_HOST_DEVICE
    attach_populator(const value_type invalid = invalid_value<value_type>())
        : kInvalid(invalid) {}

    value_type kInvalid;

    using bare_value = value_type;
    using store_value = vector_type<bare_value>;
    using serialized_storage = jagged_vector_type<bare_value>;

    using vector_view_t = vecmem::data::jagged_vector_view<bare_value>;
    using vector_buffer_t = vecmem::data::jagged_vector_buffer<bare_value>;
    using buffer_size_t = std::vector<typename vector_view_t::size_type>;

    /** Add a new value to the stored value - for host vector
     *
     * @param stored the stored value for the population
     * @param bvalue the new value to be added
     **/
    DETRAY_HOST
    void operator()(store_value &stored, bare_value &&bvalue) const {
        stored.push_back(bvalue);
        if (kSORT) {
            std::sort(stored.begin(), stored.end());
        }
    }

    /** Add a new value to the stored value - for device vector
     *
     * @param stored the stored value for the population
     * @param bvalue the new value to be added
     **/
#if defined(__CUDACC__)  // to resolve ambiguoty from host side
    DETRAY_DEVICE
    void operator()(store_value stored, bare_value &&bvalue) const {
        stored.push_back(bvalue);
    }
#endif

    /** Create a sequence of bare values, independent of the store_value.
     *
     * @param stored the stored value
     *
     * @return a sequence of bare values
     */
    DETRAY_HOST_DEVICE
    vector_type<bare_value> sequence(store_value &stored) const {
        return stored;
    }

    /** Shift operation for unified memory block
     *
     * @param stored the stored value
     * @param offset is the shift offset
     *
     **/
    DETRAY_HOST_DEVICE
    void shift(store_value &stored, const bare_value &offset) const {
        std::for_each(stored.begin(), stored.end(),
                      [&](auto &d) { d += offset; });
    }

    /** Return an initialized bin value
     **/
    DETRAY_HOST_DEVICE
    store_value init() const { return {}; }

    /** Return a vector view
     **/
    DETRAY_HOST
    static vector_view_t get_data(serialized_storage &data,
                                  vecmem::memory_resource &resource) {
        return vecmem::get_data(data, &resource);
    }
};

// convinient declaration for host replace populator
template <typename value_type = dindex>
using host_replace_populator = replace_populator<value_type>;

// convinient declaration for device replace populator
template <typename value_type = dindex>
using device_replace_populator =
    replace_populator<value_type, vecmem::device_vector>;

// convinient declaration for host complete populator
template <unsigned int kDIM, bool kSORT = false, typename value_type = dindex,
          template <typename, unsigned int> class array_type = darray>
using host_complete_populator =
    complete_populator<kDIM, kSORT, value_type, array_type>;

// convinient declaration for device complete populator
template <unsigned int kDIM, bool kSORT = false, typename value_type = dindex,
          template <typename, unsigned int> class array_type = darray>
using device_complete_populator =
    complete_populator<kDIM, kSORT, value_type, array_type,
                       vecmem::device_vector>;

// convinient declaration for host attach populator
template <bool kSORT = false, typename value_type = dindex>
using host_attach_populator = attach_populator<kSORT, value_type>;

// convinient declaration for device attach populator
template <bool kSORT = false, typename value_type = dindex>
using device_attach_populator =
    attach_populator<kSORT, value_type, vecmem::device_vector,
                     vecmem::jagged_device_vector>;

}  // namespace detray
