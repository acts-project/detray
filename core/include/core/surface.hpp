#pragma once

#include <memory>

namespace detray
{
    /** Templated surface class
     * 
     * @tparam transform_type the type of the global 3D to local 3D frame
     * @tparam source_type the type of the Source representation 
     */
    template <typename transform_type, typename source_type = int>
    class surface
    {
    public:
        /** Only allowed parameter constructor
         * 
         * @param trf the transform for positioning and 3D local frame 
         * @param src the source object this surface is representing
         * 
         **/
        surface(transform_type &&trf, source_type &&src)
            : _trf(std::move(trf)), _src(std::move(src))
        {
        }        

        ~surface() = default;
        surface(const surface& lhs) = default;
        surface() = delete;

        /** Return the transform type */
        const transform_type& transform() const 
        { return _trf; }

        /** Return the transform type */
        const source_type& source() const 
        { return _src; }

    private:
        transform_type _trf;
        source_type _src;

    };
} // namespace detray
