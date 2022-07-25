/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <memory>

#include "detray/definitions/detail/accessor.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/** Templated surface class for detector surfaces and portals
 *
 * @tparam mask_registry_t the type of the mask link representation
 * @tparam material_registry_t the type of the material link representation
 * @tparam transform_link_t the type of the transform link representation
 * @tparam volume_link_t the typ eof the volume link representation
 * @tparam source_link_t the type of the source link representation
 */
template <typename mask_regsitry_t, typename material_registry_t,
          typename transform_link_t = dindex, typename volume_link_t = dindex,
          typename source_link_t = bool>
class surface {

    public:
    // Broadcast the type of links
    using transform_link = transform_link_t;
    using mask_defs = mask_regsitry_t;
    using mask_link = typename mask_defs::link_type;
    // At least one mask type is present in any geometry
    using edge_type = typename mask_defs::template get_type<mask_defs::to_id(
        0)>::type::links_type;
    using material_defs = material_registry_t;
    using material_link = typename material_defs::link_type;
    using volume_link = volume_link_t;
    using source_link = source_link_t;

    /** Constructor with full arguments - move semantics
     *
     * @param trf the transform for positioning and 3D local frame
     * @param msk the mask/mask link for this surface
     * @param vol the volume link for this surface
     * @param src the source object/source link this surface is representing
     * @param is_pt remember whether this is a portal or not
     *
     **/
    surface(transform_link &&trf, mask_link &&mask, material_link &&material,
            volume_link &&vol, source_link &&src, bool is_pt)
        : _trf(std::move(trf)),
          _mask(std::move(mask)),
          _material(std::move(material)),
          _vol(std::move(vol)),
          _src(std::move(src)),
          _is_portal(std::move(is_pt)) {}

    /** Constructor with full arguments - copy semantics
     *
     * @param trf the transform for positioning and 3D local frame
     * @param msk the mask/mask link for this surface
     * @param vol the volume link for this surface
     * @param src the source object/source link this surface is representing
     * @param is_pt remember whether this is a portal or not
     *
     **/
    surface(const transform_link &trf, const mask_link &mask,
            const material_link &material, const volume_link vol,
            const source_link &src, bool is_pt)
        : _trf(trf),
          _mask(mask),
          _material(material),
          _vol(vol),
          _src(src),
          _is_portal(is_pt) {}

    // Portal vs module decision must be made explicitly
    surface() = default;
    surface(const surface &lhs) = default;
    ~surface() = default;

    /** Equality operator
     *
     * @param rhs is the right hand side to be compared to
     */
    DETRAY_HOST_DEVICE
    bool operator==(const surface &rhs) const {
        return (_trf == rhs._trf and _mask == rhs._mask and _vol == rhs._vol and
                _src == rhs._src and _is_portal == rhs._is_portal);
    }

    /** Update the transform index
     *
     * @param offset update the position when move into new collection
     */
    DETRAY_HOST
    void update_transform(dindex offset) { _trf += offset; }

    /** Access to the transform index */
    DETRAY_HOST_DEVICE
    const transform_link &transform() { return _trf; }

    /** @return the transform index */
    DETRAY_HOST_DEVICE
    const transform_link &transform() const { return _trf; }

    /// Update the mask link
    ///
    /// @param offset update the position when move into new collection
    DETRAY_HOST
    void update_mask(dindex offset) { _mask += offset; }

    /** Access to the mask  */
    DETRAY_HOST_DEVICE
    const mask_link &mask() { return _mask; }

    /** @return the mask link */
    DETRAY_HOST_DEVICE
    const mask_link &mask() const { return _mask; }

    /** Access to the mask id */
    DETRAY_HOST_DEVICE
    auto mask_type() { return detail::get<0>(_mask); }

    /** @return the mask link */
    DETRAY_HOST_DEVICE
    auto mask_type() const { return detail::get<0>(_mask); }

    /** Access to the mask  */
    DETRAY_HOST_DEVICE
    const auto &mask_range() { return detail::get<1>(_mask); }

    /** @return the mask link */
    DETRAY_HOST_DEVICE
    const auto &mask_range() const { return detail::get<1>(_mask); }

    /// Update the material link
    ///
    /// @param offset update the position when move into new collection
    DETRAY_HOST
    void update_material(dindex offset) { _material += offset; }

    /// Access to the material
    DETRAY_HOST_DEVICE
    const material_link &material() { return _material; }

    /// @return the material link
    DETRAY_HOST_DEVICE
    const material_link &material() const { return _material; }

    /// Access to the material id
    DETRAY_HOST_DEVICE
    auto material_type() { return detail::get<0>(_material); }

    /// @return the material link
    DETRAY_HOST_DEVICE
    auto material_type() const { return detail::get<0>(_material); }

    /// Access to the material
    DETRAY_HOST_DEVICE
    const auto &material_range() { return detail::get<1>(_material); }

    /// @return the material link
    DETRAY_HOST_DEVICE
    const auto &material_range() const { return detail::get<1>(_material); }

    /** Access to the volume */
    DETRAY_HOST_DEVICE
    volume_link volume() { return _vol; }

    /** @return the volume index */
    DETRAY_HOST_DEVICE
    const volume_link volume() const { return _vol; }

    /** @return the source link */
    DETRAY_HOST_DEVICE
    const source_link &source() const { return _src; }

    /** @return if the surface belongs to grid **/
    auto get_grid_status() const { return _in_grid; }

    /** set if the surface belongs to grid **/
    void set_grid_status(bool status) { _in_grid = status; }

    /** Is this instance a portal in the sense of the unified_index_geometry? */
    DETRAY_HOST_DEVICE
    bool is_portal() const { return _is_portal; }

    private:
    transform_link_t _trf;
    mask_link _mask;
    material_link _material;
    volume_link_t _vol;
    source_link_t _src;
    bool _is_portal;
    bool _in_grid = false;
};

}  // namespace detray
