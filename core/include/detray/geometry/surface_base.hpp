/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <memory>

#include "detray/utils/indexing.hpp"

namespace detray {
/** Templated surface base class for detector surfaces and portals
 *
 * @tparam transform_link the type of the transform/transform link for global 3D
 * to local 3D frame
 * @tparam mask_link the type of the mask/mask link representation
 * @tparam volume_link the typ eof the volume/volume link representation
 * @tparam source_link the type of the source/source link representation
 */
template <typename transform_link, typename mask_link = dindex,
          typename volume_link = dindex, typename source_link = bool,
          typename edge_link = dindex>
class surface_base {
    public:
    // Broadcast the type of links
    using transform_links = transform_link;
    using mask_links = mask_link;
    using volume_links = volume_link;
    using source_links = source_link;
    using edge_links = edge_link;

    /** Constructor with full arguments - move semantics
     *
     * @param trf the transform for positioning and 3D local frame
     * @param msk the mask/mask link for this surface
     * @param vol the volume link for this surface
     * @param src the source object/source link this surface is representing
     *
     **/
    surface_base(transform_link &&trf, mask_link &&mask, volume_link &&vol,
                 source_link &&src)
        : _trf(std::move(trf)),
          _mask(std::move(mask)),
          _vol(std::move(vol)),
          _src(std::move(src)) {}

    /** Constructor with full arguments - copy semantics
     *
     * @param trf the transform for positioning and 3D local frame
     * @param msk the mask/mask link for this surface
     * @param vol the volume link for this surface
     * @param src the source object/source link this surface is representing
     *
     **/
    surface_base(const transform_link &trf, const mask_link &mask,
                 const volume_link &vol, const source_link &src)
        : _trf(trf), _mask(mask), _vol(vol), _src(src) {}

    ~surface_base() = default;
    surface_base(const surface_base &lhs) = default;
    surface_base() = delete;

    /** Equality operator
     *
     * @param rhs is the right hand side to be compared to
     */
    bool operator==(
        const surface_base<transform_link, mask_link, source_link> &rhs) const {
        return (_trf == rhs.__trf and _mask == rhs._mask and
                _vol == rhs._vol and _src == rhs._src and _edg == rhs._edg);
    }

    /** Explicitly set edge, since not all geometries keep the links here */
    void set_edge(const edge_link &edg) { _edg = edg; }

    /** Access to the edge information (next volume etc.)  */
    const edge_link &edge() const { return _edg; }

    /** Return the edge information (next volume etc.)  */
    edge_link &edge() { return _edg; }

    /** Return the transform type */
    transform_link &transform() { return _trf; }

    /** Return the transform type */
    const transform_link &transform() const { return _trf; }

    /** Access to the mask  */
    mask_link &mask() { return _mask; }

    /** Return the mask */
    const mask_link &mask() const { return _mask; }

    /** Access to the volume */
    volume_link &volume() { return _vol; }

    /** Return the mask */
    const volume_link &volume() const { return _vol; }

    /** Return the source/source link type */
    const source_link &source() const { return _src; }

    private:
    transform_link _trf;
    mask_link _mask;
    volume_link _vol;
    source_link _src;
    edge_link _edg = {};
};

}  // namespace detray
