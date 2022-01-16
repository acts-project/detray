/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <memory>

#include "detray/definitions/qualifiers.hpp"
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
    surface_base() = default;

    /** Equality operator
     *
     * @param rhs is the right hand side to be compared to
     */
    DETRAY_HOST_DEVICE
    bool operator==(const surface_base &rhs) const {
        return (_trf == rhs._trf and _mask == rhs._mask and _vol == rhs._vol and
                _src == rhs._src and _edg == rhs._edg);
    }

    /** Explicitly set edge, since not all geometries keep the links here */
    DETRAY_HOST_DEVICE
    void set_edge(const edge_link &edg) { _edg = edg; }

    /** Access to the edge information (next volume etc.)  */
    DETRAY_HOST_DEVICE
    const edge_link &edge() const { return _edg; }

    /** @return the edge information (next volume etc.)  */
    DETRAY_HOST_DEVICE
    edge_link &edge() { return _edg; }

    /** Access to the transform index */
    DETRAY_HOST_DEVICE
    transform_link &transform() { return _trf; }

    /** @return the transform index */
    DETRAY_HOST_DEVICE
    const transform_link &transform() const { return _trf; }

    /** Access to the mask  */
    DETRAY_HOST_DEVICE
    mask_link &mask() { return _mask; }

    /** @return the mask link */
    DETRAY_HOST_DEVICE
    const mask_link &mask() const { return _mask; }

    /** Access to the volume */
    DETRAY_HOST_DEVICE
    volume_link &volume() { return _vol; }

    /** @return the volume index */
    DETRAY_HOST_DEVICE
    const volume_link &volume() const { return _vol; }

    /** @return the source link */
    DETRAY_HOST_DEVICE
    const source_link &source() const { return _src; }

    /** get if the surface belongs to grid **/
    auto get_grid_status() const { return _in_grid; }

    /** set if the surface belongs to grid **/
    void set_grid_status(bool status) { _in_grid = status; }

    /** Is this instance a portal in the sense of the unified_index_geometry? */
    DETRAY_HOST_DEVICE
    bool is_portal() const { return not(std::get<0>(_edg) == _vol); }

    private:
    transform_link _trf;
    mask_link _mask;
    volume_link _vol;
    source_link _src;
    edge_link _edg = {};
    bool _in_grid = false;
};

}  // namespace detray
