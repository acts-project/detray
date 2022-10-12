
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

namespace {

enum surface_id { e_sensitive, e_portal, e_passive };
}

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
          typename transform_link_t = dindex, typename source_link_t = bool>
class surface {

    public:
    /// Link type of the mask to a volume.
    using volume_link_type = dindex;
    // Broadcast the type of links
    using transform_link = transform_link_t;
    /// might be a single mask, a range of masks or a multiindex in the future
    using mask_link = mask_link_t;
    using mask_id = typename mask_link::id_type;
    using material_link = material_link_t;
    using material_id = typename material_link::id_type;
    ;
    using source_link = source_link_t;

    /// Constructor with full arguments - move semantics
    ///
    /// @param trf the transform for positioning and 3D local frame
    /// @param mask the type and index of the mask for this surface
    /// @param material the type and index of the material for this surface
    /// @param vol the volume this surface belongs to
    /// @param src the source object/source link this surface is representing
    /// @param is_portal remember whether this is a portal or not
    surface(transform_link &&trf, mask_link &&mask, material_link &&material,
            dindex volume, source_link &&src, bool is_portal)
        : _trf(std::move(trf)),
          _mask(std::move(mask)),
          _material(std::move(material)),
          _volume(volume),
          _src(std::move(src)),
          _is_portal(is_portal) {}

    /// Constructor with full arguments - copy semantics
    ///
    /// @param trf the transform for positioning and 3D local frame
    /// @param mask the type and index of the mask for this surface
    /// @param material the type and index of the material for this surface
    /// @param vol the volume this surface belongs to
    /// @param src the source object/source link this surface is representing
    /// @param is_portal remember whether this is a portal or not
    surface(const transform_link trf, const mask_link &mask,
            const material_link &material, const dindex volume,
            const source_link &src, bool is_portal)
        : _trf(trf),
          _mask(mask),
          _material(material),
          _volume(volume),
          _src(src),
          _is_portal(is_portal) {}

    /// Portal vs module decision must be made explicitly
    surface() = default;
    surface(const surface &lhs) = default;
    ~surface() = default;

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE
    auto operator==(const surface &rhs) const -> bool {
        return (_trf == rhs._trf and _mask == rhs._mask and
                _volume == rhs._volume and _src == rhs._src and
                _is_portal == rhs._is_portal);
    }

    /// Update the transform index
    ///
    /// @param offset update the position when move into new collection
    DETRAY_HOST
    auto update_transform(dindex offset) -> void { _trf += offset; }

    /// Access to the transform index
    DETRAY_HOST_DEVICE
    auto transform() -> const transform_link & { return _trf; }

    /// @return the transform index
    DETRAY_HOST_DEVICE
    auto transform() const -> const transform_link & { return _trf; }

    /// Update the mask link
    ///
    /// @param offset update the position when move into new collection
    DETRAY_HOST
    auto update_mask(dindex offset) -> void { _mask += offset; }

    /// Access to the mask
    DETRAY_HOST_DEVICE
    auto mask() -> const mask_link & { return _mask; }

    /// @return the mask link
    DETRAY_HOST_DEVICE
    auto mask() const -> const mask_link & { return _mask; }

    /// Access to the mask id
    DETRAY_HOST_DEVICE
    auto mask_type() -> typename mask_link::id_type {
        return detail::get<0>(_mask);
    }

    /// @return the mask link
    auto mask_type() const -> typename mask_link::id_type {
        return detail::get<0>(_mask);
    }

    /// Access to the mask
    DETRAY_HOST_DEVICE
    auto mask_range() -> typename mask_link::index_type & {
        return detail::get<1>(_mask);
    }

    /// @return the mask link
    DETRAY_HOST_DEVICE
    auto mask_range() const -> const typename mask_link::index_type & {
        return detail::get<1>(_mask);
    }

    /// Update the material link
    ///
    /// @param offset update the position when move into new collection
    DETRAY_HOST
    auto update_material(dindex offset) -> void { _material += offset; }

    /// Access to the material
    DETRAY_HOST_DEVICE
    auto material() -> const material_link & { return _material; }

    /// @return the material link
    DETRAY_HOST_DEVICE
    auto material() const -> const material_link & { return _material; }

    /// Access to the material id
    DETRAY_HOST_DEVICE
    auto material_type() -> typename material_link::id_type & {
        return detail::get<0>(_material);
    }

    /// @return the material link
    DETRAY_HOST_DEVICE
    auto material_type() const -> const typename material_link::id_type & {
        return detail::get<0>(_material);
    }

    /// Access to the material
    DETRAY_HOST_DEVICE
    auto material_range() -> typename material_link::index_type & {
        return detail::get<1>(_material);
    }

    /// @return the material link
    DETRAY_HOST_DEVICE
    auto material_range() const -> const typename material_link::index_type & {
        return detail::get<1>(_material);
    }

    /// Access to the volume
    DETRAY_HOST_DEVICE
    auto volume() -> dindex { return _volume; }

    /// @return the volume index
    DETRAY_HOST_DEVICE
    auto volume() const -> dindex { return _volume; }

    /// @return the source link
    DETRAY_HOST_DEVICE
    auto source() const -> const source_link & { return _src; }

    /// Is this instance a portal in the sense of the unified_index_geometry?
    DETRAY_HOST_DEVICE
    auto is_portal() const -> bool { return _is_portal; }

    private:
    transform_link_t _trf{};
    mask_link _mask{};
    material_link _material{};
    dindex _volume{dindex_invalid};
    source_link_t _src{};
    bool _is_portal = false;
};

}  // namespace detray