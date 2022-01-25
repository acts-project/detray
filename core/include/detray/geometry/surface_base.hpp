/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <memory>

#include "detray/definitions/detail/accessor.hpp"
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
template <typename mask_regsitry_t, typename intersection_kernel_t,
          typename transform_link_t = dindex, typename volume_link_t = dindex,
          typename source_link_t = bool,
          typename edge_link_t = typed_index<unsigned int, dindex>>
class surface_base {
    public:
    // Broadcast the type of links
    using transform_link = transform_link_t;
    using mask_defs = mask_regsitry_t;
    using mask_link = typename mask_defs::link_type;
    using volume_link = volume_link_t;
    using source_link = source_link_t;
    using edge_link = edge_link_t;

    /** Constructor with full arguments - move semantics
     *
     * @param trf the transform for positioning and 3D local frame
     * @param msk the mask/mask link for this surface
     * @param vol the volume link for this surface
     * @param src the source object/source link this surface is representing
     *
     **/
    surface_base(transform_link_t &&trf, mask_link &&mask, volume_link_t &&vol,
                 source_link_t &&src)
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
    surface_base(const transform_link_t &trf, const mask_link &mask,
                 const volume_link_t &vol, const source_link_t &src)
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
    void set_edge(const edge_link_t &edg) { _edg = edg; }

    /** Access to the edge information (next volume etc.)  */
    DETRAY_HOST_DEVICE
    const edge_link_t &edge() const { return _edg; }

    /** @return the edge information (next volume etc.)  */
    DETRAY_HOST_DEVICE
    edge_link_t &edge() { return _edg; }

    /** Access to the transform index */
    DETRAY_HOST_DEVICE
    transform_link_t &transform() { return _trf; }

    /** @return the transform index */
    DETRAY_HOST_DEVICE
    const transform_link_t &transform() const { return _trf; }

    /** Access to the mask  */
    DETRAY_HOST_DEVICE
    mask_link &mask() { return _mask; }

    /** @return the mask link */
    DETRAY_HOST_DEVICE
    const mask_link &mask() const { return _mask; }

    /** Access to the mask  */
    DETRAY_HOST_DEVICE
    auto &mask_type() { return detail::get<0>(_mask); }

    /** @return the mask link */
    DETRAY_HOST_DEVICE
    const auto &mask_type() const { return detail::get<0>(_mask); }

    /** Access to the mask  */
    DETRAY_HOST_DEVICE
    auto &mask_range() { return detail::get<1>(_mask); }

    /** @return the mask link */
    DETRAY_HOST_DEVICE
    const auto &mask_range() const { return detail::get<1>(_mask); }

    /** Access to the volume */
    DETRAY_HOST_DEVICE
    volume_link_t &volume() { return _vol; }

    /** @return the volume index */
    DETRAY_HOST_DEVICE
    const volume_link_t &volume() const { return _vol; }

    /** @return the source link */
    DETRAY_HOST_DEVICE
    const source_link_t &source() const { return _src; }

    /** get if the surface belongs to grid **/
    auto get_grid_status() const { return _in_grid; }

    /** set if the surface belongs to grid **/
    void set_grid_status(bool status) { _in_grid = status; }

    /** Is this instance a portal in the sense of the unified_index_geometry? */
    DETRAY_HOST_DEVICE
    bool is_portal() const { return not(detail::get<0>(_edg) == _vol); }

    /** Kernel method that updates the intersections
     *
     * @tparam track_t The type of the track/context
     * @tparam transform_container_t The type of the transform container
     * @tparam mask_container_t The type of the mask container
     *
     * @param track the track information including the contexts
     * @param contextual_transform the transform container
     * @param masks the tuple mask container to for the intersection
     *
     * @return  an intersection struct (invalid if no intersection was found)
     **/
    template <typename track_t, typename transform_container_t,
              typename mask_container_t>
    inline const auto intersect(
        const track_t &track,
        const transform_container_t &contextual_transforms,
        const mask_container_t &masks) const {
        return intersection_kernel_t{}(*this, track, contextual_transforms,
                                       masks);
    }

    private:
    transform_link_t _trf;
    mask_link _mask;
    volume_link_t _vol;
    source_link_t _src;
    edge_link_t _edg = {};
    bool _in_grid = false;
};

}  // namespace detray
