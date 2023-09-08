/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/barcode.hpp"

// Sysytem include(s)
#include <memory>

namespace detray {

/// Templated surface class for detector surfaces and portals.
///
/// @note might be holding multiple surfaces in the future
///
/// @tparam mask_regsitry_t the type collection of masks that can be linked
///                         to the surface
/// @tparam material_registry_t the type collection of material that can be
///                             linked to the surface
/// @tparam transform_link_t how to reference the surfaces transforms
/// @tparam source_link_t the type of the source link representation
template <typename mask_link_t = dtyped_index<dindex, dindex>,
          typename material_link_t = dtyped_index<dindex, dindex>,
          typename transform_link_t = dindex,
          typename navigation_link_t = dindex, typename source_link_t = dindex>
class surface_descriptor {

    public:
    /// Link type of the mask to a volume.
    using navigation_link = navigation_link_t;
    // Broadcast the type of links
    using transform_link = transform_link_t;
    /// might be a single mask, a range of masks or a multiindex in the future
    using mask_link = mask_link_t;
    using mask_id = typename mask_link::id_type;
    using material_link = material_link_t;
    using material_id = typename material_link::id_type;
    using source_link = source_link_t;

    /// Constructor with full arguments - move semantics
    ///
    /// @param trf the transform for positioning and 3D local frame
    /// @param mask the type and index of the mask for this surface
    /// @param material the type and index of the material for this surface
    /// @param vol the volume this surface belongs to
    /// @param src the source object/source link this surface is representing
    /// @param sf_id remember whether this is a portal or not
    DETRAY_HOST
    constexpr surface_descriptor(transform_link &&trf, mask_link &&mask,
                                 material_link &&material, dindex volume,
                                 source_link &&src, surface_id sf_id)
        : _mask(std::move(mask)),
          _material(std::move(material)),
          _trf(std::move(trf)),
          _src(std::move(src)) {

        m_barcode = geometry::barcode{}.set_volume(volume).set_id(sf_id);
    }

    /// Constructor with full arguments - copy semantics
    ///
    /// @param trf the transform for positioning and 3D local frame
    /// @param mask the type and index of the mask for this surface
    /// @param material the type and index of the material for this surface
    /// @param vol the volume this surface belongs to
    /// @param src the source object/source link this surface is representing
    /// @param sf_id remember whether this is a portal or not
    DETRAY_HOST
    constexpr surface_descriptor(const transform_link trf,
                                 const mask_link &mask,
                                 const material_link &material,
                                 const dindex volume, const source_link &src,
                                 surface_id sf_id)
        : _mask(mask), _material(material), _trf(trf), _src(src) {
        m_barcode = geometry::barcode{}.set_volume(volume).set_id(sf_id);
    }

    constexpr surface_descriptor() = default;
    surface_descriptor(const surface_descriptor &lhs) = default;
    ~surface_descriptor() = default;

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const surface_descriptor &rhs) const -> bool {
        return (_mask == rhs._mask and _material == rhs._material and
                _trf == rhs._trf and _src == rhs._src and
                m_barcode == rhs.m_barcode);
    }

    /// Sets a new surface barcode
    DETRAY_HOST_DEVICE
    auto set_barcode(const geometry::barcode bcd) -> void { m_barcode = bcd; }

    /// @returns the surface barcode
    DETRAY_HOST_DEVICE
    constexpr auto barcode() const -> geometry::barcode { return m_barcode; }

    /// Sets a new surface id (portal/passive/sensitive)
    DETRAY_HOST_DEVICE
    auto set_id(const surface_id new_id) -> void { m_barcode.set_id(new_id); }

    /// @returns the surface id (sensitive, passive or portal)
    DETRAY_HOST_DEVICE
    constexpr auto id() const -> surface_id { return m_barcode.id(); }

    /// Sets a new volume link (index in volume collection of detector)
    DETRAY_HOST
    auto set_volume(const dindex new_idx) -> void {
        m_barcode.set_volume(new_idx);
    }

    /// @returns the surface id (sensitive, passive or portal)
    DETRAY_HOST_DEVICE
    constexpr auto volume() const -> dindex { return m_barcode.volume(); }

    /// Sets a new surface index (index in surface collection of surface store)
    DETRAY_HOST_DEVICE
    auto set_index(const dindex new_idx) -> void {
        m_barcode.set_index(new_idx);
    }

    /// @returns the surface id (sensitive, passive or portal)
    DETRAY_HOST_DEVICE
    constexpr auto index() const -> dindex { return m_barcode.index(); }

    /// Update the transform index
    ///
    /// @param offset update the position when move into new collection
    DETRAY_HOST
    auto update_transform(dindex offset) -> void { _trf += offset; }

    /// @return the transform index
    DETRAY_HOST_DEVICE
    constexpr auto transform() const -> const transform_link & { return _trf; }

    /// Update the mask link
    ///
    /// @param offset update the position when move into new collection
    DETRAY_HOST
    auto update_mask(dindex offset) -> void { _mask += offset; }

    /// @return the mask link
    DETRAY_HOST_DEVICE
    constexpr auto mask() const -> const mask_link & { return _mask; }

    /// Update the material link
    ///
    /// @param offset update the position when move into new collection
    DETRAY_HOST
    auto update_material(dindex offset) -> void { _material += offset; }

    /// Access to the material
    DETRAY_HOST_DEVICE
    constexpr auto material() -> material_link & { return _material; }

    /// @return the material link
    DETRAY_HOST_DEVICE
    constexpr auto material() const -> const material_link & {
        return _material;
    }

    /// @return the source link
    DETRAY_HOST_DEVICE
    constexpr auto source() const -> const source_link & { return _src; }

    /// @returns true if the surface is a senstive detector module.
    DETRAY_HOST_DEVICE
    constexpr auto is_sensitive() const -> bool {
        return m_barcode.id() == surface_id::e_sensitive;
    }

    /// @returns true if the surface is a portal.
    DETRAY_HOST_DEVICE
    constexpr auto is_portal() const -> bool {
        return m_barcode.id() == surface_id::e_portal;
    }

    /// @returns true if the surface is a passive detector element.
    DETRAY_HOST_DEVICE
    constexpr auto is_passive() const -> bool {
        return m_barcode.id() == surface_id::e_passive;
    }

    private:
    geometry::barcode m_barcode;
    mask_link _mask;
    material_link _material;
    transform_link_t _trf;
    source_link_t _src;
};

}  // namespace detray
