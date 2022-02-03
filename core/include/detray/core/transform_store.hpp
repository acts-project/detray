/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <vecmem/memory/memory_resource.hpp>

#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/enumerate.hpp"
#include "detray/utils/indexing.hpp"

namespace detray {

using transform3 = __plugin::transform3<detray::scalar>;

/** A static inplementation of an alignable transform store */
template <template <typename...> class vector_t = dvector,
          typename context_t = dindex>
class static_transform_store {
    public:
    using storage = vector_t<transform3>;
    using context = context_t;

    static_transform_store() = default;

    /** Constructor with vecmem memory resource
     **/
    static_transform_store(vecmem::memory_resource &resource)
        : _data(&resource) {}

    /** Constructor from static_transform_store_data
     **/
    template <typename transform_store_data_t,
              std::enable_if_t<!std::is_base_of_v<vecmem::memory_resource,
                                                  transform_store_data_t>,
                               bool> = true>
    DETRAY_DEVICE static_transform_store(transform_store_data_t &store_data)
        : _data(store_data._data) {}

    /** Elementwise access. Needs []operator for storage type for now */
    DETRAY_HOST_DEVICE
    inline auto operator[](const dindex i) { return _data[i]; }
    DETRAY_HOST_DEVICE
    inline auto operator[](const dindex i) const { return _data[i]; }

    /** Forward iterator : Contextual STL like API
     *
     * @param ctx The context of the call (ignored)
     */
    DETRAY_HOST_DEVICE
    auto begin(const context & /*ctx*/) const { return _data.begin(); }

    /** Forward iterator : Contextual STL like API
     *
     * @param ctx The context of the call (ignored)
     */
    DETRAY_HOST_DEVICE
    auto end(const context & /*ctx*/) const { return _data.end(); }

    /** Access to a predefined range of elements
     *
     * @param start start index of rage
     * @param end end index of rage
     * @param ctx The context of the call (ignored)
     *
     * @return range restricted iterator
     */
    DETRAY_HOST_DEVICE
    const inline auto range(const size_t begin, const size_t end,
                            const context & /*ctx*/) const {
        return iterator_range(_data, dindex_range{begin, end});
    }

    /** Access to a predefined range of elements
     *
     * @param range index range in data store
     * @param ctx The context of the call (ignored)
     *
     * @return range restricted iterator
     */
    template <typename range_t>
    DETRAY_HOST_DEVICE const inline auto range(range_t &&range,
                                               const context & /*ctx*/) const {
        return iterator_range(_data, std::forward<range_t>(range));
    }

    /** Reserve memory : Contextual STL like API
     *
     * @param ctx The context of the call (ignored)
     * @param tidx The size of the reserved container memory
     */
    DETRAY_HOST
    void reserve(const context & /*ctx*/, size_t n_size) {
        _data.reserve(n_size);
    }

    /** Emplace back : Contextual STL like API
     *
     * @param ctx The context of the call (ignored)
     * @param args Constructor arguments
     */
    template <class... Args>
    DETRAY_HOST auto &emplace_back(const context & /*ctx*/, Args &&... args) {
        return _data.emplace_back(std::forward<Args>(args)...);
    }

    /** Push back : Contextual STL like API, copy semantics
     *
     * @param ctx The context of the call (ignored)
     * @param tcf The transform to be filled
     */
    DETRAY_HOST
    void push_back(const context & /*ctx*/, const transform3 &tf) {
        _data.push_back(tf);
    }

    /** Push back : Contextual STL like API, move semantics
     *
     * @param ctx The context of the call (ignored)
     * @param tcf The transform to be filled
     */
    DETRAY_HOST
    void push_back(const context & /*ctx*/, const transform3 &&tf) {
        _data.push_back(std::move(tf));
    }

    /** Size : Contextual STL like API
     * @param ctx The context of the call (ignored)
     */
    DETRAY_HOST_DEVICE
    size_t size(const context & /*ctx*/) const { return _data.size(); }

    /** Empty : Contextual STL like API
     * @param ctx The context of the call (ignored)
     */
    DETRAY_HOST_DEVICE
    bool empty(const context & /*ctx*/) { return _data.empty(); }

    /** Empty : Contextual STL like API
     * @param ctx The context of the call (ignored)
     */
    DETRAY_HOST_DEVICE
    bool empty(const context & /*ctx*/) const { return _data.empty(); }

    /** Append a transform store to an existing one
     *
     * @param ctx The context of the call (ignored)
     * @param other The transform store, move semantics
     *
     * @note in general can throw an exception
     */
    DETRAY_HOST
    void append(const context &ctx,
                static_transform_store<vector_t> &&other) noexcept(false) {
        _data.reserve(_data.size() + other.size(ctx));
        _data.insert(_data.end(), std::make_move_iterator(other.begin(ctx)),
                     std::make_move_iterator(other.end(ctx)));
    }

    /** Add a new bunch of transforms for a new context - move semantics
     *
     * @param ctx The context of the call (ignored)
     * @param trfs The transform container, move semantics
     *
     * @note in general can throw an exception
     */
    DETRAY_HOST
    void add_contextual_transforms(const context & /*ctx*/,
                                   storage &&trfs) noexcept(false) {
        _data = std::move(trfs);
    }

    /** Append a bunch of transforms to existing context - move semantics
     *
     * @param ctx The context of the call (ignored)
     * @param trfs The transform container, move semantics
     *
     * @note in general can throw an exception
     */
    DETRAY_HOST
    void append_contextual_transforms(const context & /*ctx*/,
                                      storage &&trfs) noexcept(false) {
        _data.reserve(_data.size() + trfs.size());
        _data.insert(_data.end(), std::make_move_iterator(trfs.begin()),
                     std::make_move_iterator(trfs.end()));
    }

    /** Get the contextual transform
     *
     * @param ctx The context of the call (ignored)
     * @param tidx The transform index
     */
    DETRAY_HOST_DEVICE
    const transform3 &contextual_transform(const context & /*ctx*/,
                                           dindex tidx) const {
        return _data[tidx];
    }

    DETRAY_HOST_DEVICE
    vector_t<transform3> &data() { return _data; }

    private:
    /** Common to surfaces & portals: transform store */
    vector_t<transform3> _data;
};

/** A static inplementation of transform store data for device*/
template <typename transform_store_t>
struct static_transform_store_data {

    /** Constructor from transform store
     *
     * @param store is the input transform store data from host
     **/
    static_transform_store_data(transform_store_t &store)
        : _data(vecmem::get_data(store.data())) {}

    vecmem::data::vector_view<transform3> _data;
};

/** Get transform_store_data
 **/
template <template <typename...> class vector_t>
inline static_transform_store_data<static_transform_store<vector_t> > get_data(
    static_transform_store<vector_t> &store) {
    return static_transform_store_data(store);
}

}  // namespace detray
