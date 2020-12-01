/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <memory>

namespace detray
{
    /** Templated surface class
     * 
     * @tparam transform_type the type of the global 3D to local 3D frame
     * @tparam mask_link the type of the maks/maskj link representation
     * @tparam source_link the type of the source/source link representation 
     */
    template <typename transform_type, typename mask_link = int, typename source_link = int>
    class surface
    {        
    public:
        /** Broadcast the transform type */
        using transform3 = transform_type;

        /** Constructor with full arguments
         * 
         * @param trf the transform for positioning and 3D local frame 
         * @param msk thie mask/mask link for this surface 
         * @param src the source object/source link this surface is representing
         * 
         **/
        surface(transform_type &&trf, mask_link&& mask, source_link &&src)
            : _trf(std::move(trf)), _mask(mask), _src(std::move(src))
        {
        }        

        ~surface() = default;
        surface(const surface& lhs) = default;
        surface() = delete;

        /** Return the transform type */
        const transform_type& transform() const 
        { return _trf; }

        /** Return the transform type */
        const mask_link& mask() const 
        { return _mask; }

        /** Return the source/source link type */
        const source_link& source() const 
        { return _src; }

    private:
        transform_type _trf;
        mask_link _mask;
        source_link _src;

    };
} // namespace detray
