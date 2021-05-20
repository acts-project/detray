/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "utils/indexing.hpp"
#include "utils/containers.hpp"

namespace detray
{

    using transform3 = __plugin::transform3;

    /** A static inplementation of an alignable transform store */
    class static_transform_store
    {
    public:
        /** Empty context type struct */
        struct context
        {
        };

        using storage = dvector<transform3>;

        /** Reserve memory : Contextual STL like API
         *
         * @param ctx The context of the call (ignored)
         * @param tidx The size of the reserved container memory
         */
        void reserve(const context & /*ctx*/, size_t n_size)
        {
            _data.reserve(n_size);
        }

        /** Push back : Contextual STL like API, copy semantics
         *
         * @param ctx The context of the call (ignored)
         * @param tcf The transform to be filled
         */
        void push_back(const context & /*ctx*/, const transform3 &tf)
        {
            _data.push_back(tf);
        }

        /** Push back : Contextual STL like API, move semantics
         *
         * @param ctx The context of the call (ignored)
         * @param tcf The transform to be filled
         */
        void push_back(const context & /*ctx*/, const transform3 &&tf)
        {
            _data.push_back(std::move(tf));
        }

        /** Size : Contextual STL like API
         * @param ctx The context of the call (ignored)
         */ 
        size_t size(const context & /*ctx*/)
        {
            return _data.size();
        }

        /** Empty : Contextual STL like API
         * @param ctx The context of the call (ignored)
         */
        bool empty(const context & /*ctx*/)
        {
            return _data.empty();
        }

        /** Add a new bunch of (contextual) transforms - move semantics
         *
         * @param ctx The context of the call (ignored)
         * @param trfs The transform container, move semantics
         *
         * @note in general can throw an exception
         */
        void add_contextual_transforms(const context & /*ctx*/, storage &&trfs) noexcept(false)
        {
            _data = std::move(trfs);
        }

        /** Get the contextual transform
         *
         * @param ctx The context of the call (ignored)
         * @param tidx The transform index
         */
        const transform3 &contextual_transform(const context & /*ctx*/, dindex tidx) const
        {
            return _data[tidx];
        }

    private:
        /** Common to surfaces & portals: transform store */
        dvector<transform3> _data;
    };

} // namespace detray
