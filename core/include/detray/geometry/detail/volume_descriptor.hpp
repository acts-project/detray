/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/predefined_materials.hpp"

namespace detray {

/// @brief The detray detector volume descriptor.
///
/// Contains the data and links to describe a detector volume. This is the type
/// that is stored in the detector data stores.
///
/// @tparam ID enum of object types contained in the volume
///         (@see @c detector_metadata ).
/// @tparam link_t the type of link to the volumes surfaces finder(s)
///         (accelerator structure, e.g. a grid). The surface finder types
///         cannot be given directly, since the containers differ between host
///         and device. The surface finders reside in an 'unrollable tuple
///         container' and are called per volume in the navigator during local
///         navigation.
template <typename ID, typename scalar_t,
          typename link_t = dtyped_index<dindex, dindex>>
class volume_descriptor {

    public:
    /// Ids of objects that can be distinguished by the volume
    using object_id = ID;

    /// How to access the surface ranges in the detector surface lookup
    using sf_link_type =
        dmulti_index<dindex_range, static_cast<std::size_t>(surface_id::e_all)>;

    /// How to access objects (e.g. sensitives/passives/portals) in this
    /// volume. Keeps one accelerator structure link per object type (by ID):
    ///
    /// link_t : id and index of the accelerator structure in the detector's
    ///          surface store.
    ///
    /// E.g. a 'portal' can be found under @c ID::e_portal in this link,
    /// and will then receive link to the @c brute_force_finder that holds the
    /// portals (the accelerator structure's id and index).
    using accel_link_type = dmulti_index<link_t, ID::e_size>;

    /// Default constructor builds an ~infinitely long cylinder
    constexpr volume_descriptor() = default;

    /// Constructor from shape id.
    ///
    /// @param id id values that determines how to interpret the bounds.
    explicit constexpr volume_descriptor(const volume_id id) : m_id{id} {}

    /// @returns the volume shape id, e.g. 'cylinder'
    DETRAY_HOST_DEVICE
    constexpr auto id() const -> volume_id { return m_id; }

    /// @returns the index of the volume in the detector volume container.
    DETRAY_HOST_DEVICE
    constexpr auto index() const -> dindex { return m_index; }

    /// @param index the index of the volume in the detector volume container.
    DETRAY_HOST
    constexpr auto set_index(const dindex index) -> void { m_index = index; }

    /// @returns the index of the volume tranform in the transform store.
    DETRAY_HOST_DEVICE
    constexpr auto transform() const -> dindex { return m_transform; }

    /// @param index the index of the volume in the detector volume container.
    DETRAY_HOST
    constexpr auto set_transform(const dindex trf_idx) -> void {
        m_transform = trf_idx;
    }

    /// @returns the volume material
    DETRAY_HOST_DEVICE
    constexpr auto material() const -> const detray::material<scalar_t>& {
        return m_vol_mat;
    }

    DETRAY_HOST
    constexpr auto set_material(const detray::material<scalar_t>& mat) -> void {
        m_vol_mat = mat;
    }

    /// @returns surface link for all object types - const
    DETRAY_HOST_DEVICE constexpr auto sf_link() const -> const sf_link_type& {
        return m_sf_links;
    }

    /// @returns surface descriptor link for a specific type of object - const
    template <surface_id id>
    DETRAY_HOST_DEVICE constexpr auto sf_link() const -> const
        typename sf_link_type::index_type& {
        return detail::get<static_cast<uint>(id)>(m_sf_links);
    }

    /// @returns surface descriptor link for a specific type of object
    template <surface_id id>
    DETRAY_HOST_DEVICE constexpr auto sf_link() ->
        typename sf_link_type::index_type& {
        return detail::get<static_cast<uint>(id)>(m_sf_links);
    }

    /// @returns surface descriptor link for all surface types
    DETRAY_HOST_DEVICE constexpr auto full_sf_range() const ->
        typename sf_link_type::index_type {

        using idx_range_t = typename sf_link_type::index_type;

        const auto& pt_range = sf_link<surface_id::e_portal>();
        const auto& sen_range = sf_link<surface_id::e_sensitive>();
        const auto& psv_range = sf_link<surface_id::e_passive>();

        // Portal range is never empty
        dindex min{detail::get<0>(pt_range)};
        dindex max{detail::get<1>(pt_range)};

        constexpr idx_range_t empty{};
        if (sen_range != empty) {
            min = detail::get<0>(sen_range) < min ? detail::get<0>(sen_range)
                                                  : min;
            max = detail::get<1>(sen_range) > max ? detail::get<1>(sen_range)
                                                  : max;
        }

        if (psv_range != empty) {
            min = detail::get<0>(psv_range) < min ? detail::get<0>(psv_range)
                                                  : min;
            max = detail::get<1>(psv_range) > max ? detail::get<1>(psv_range)
                                                  : max;
        }

        return idx_range_t{min, max};
    }

    /// @returns a surface index with respect to the volume surface range
    DETRAY_HOST_DEVICE constexpr dindex to_local_sf_index(dindex sf_idx) const {

        auto full_range = full_sf_range();

        assert(full_range[0] <= sf_idx);
        assert(sf_idx < full_range[1]);

        return sf_idx - full_range[0];
    }

    /// @returns a surface index with respect to the global detector containers
    DETRAY_HOST_DEVICE constexpr dindex to_global_sf_index(
        dindex sf_idx) const {

        auto full_range = full_sf_range();
        dindex glob_index{sf_idx + full_range[0]};

        assert(full_range[0] <= glob_index);
        assert(glob_index < full_range[1]);

        return glob_index;
    }

    /// Set or update the index into a geometry container identified by the
    /// obj_id.
    ///
    /// @note There is no check of overlapping index ranges between the object
    /// types. Use with care!
    ///
    /// @param other Surface index range
    template <surface_id id>
    DETRAY_HOST auto update_sf_link(
        const typename sf_link_type::index_type& other) noexcept -> void {
        auto& rg = sf_link<id>();
        // Range not set yet - initialize
        constexpr typename sf_link_type::index_type empty{};
        if (rg == empty) {
            rg = other;
        } else {
            // Update upper border
            assert(detail::get<1>(rg) == detail::get<0>(other));
            detail::get<1>(rg) = detail::get<1>(other);
        }
    }

    /// Set or update the index into a geometry container identified by the
    /// obj_id.
    ///
    /// @note There is no check of overlapping index ranges between the object
    /// types. Use with care!
    ///
    /// @param shift shift of the surface range in a larger container.
    /// @param n_surfaces the number of surfaces in this range.
    template <surface_id id>
    DETRAY_HOST auto update_sf_link(std::size_t shift,
                                    std::size_t n_surfaces = 0) noexcept
        -> void {
        auto& rg = sf_link<id>();
        // Range not set yet - initialize
        constexpr typename sf_link_type::index_type empty{};
        if (rg == empty) {
            rg = {0u, static_cast<dindex>(n_surfaces)};
        }
        // Update
        detail::get<0>(rg) += static_cast<dindex>(shift);
        detail::get<1>(rg) += static_cast<dindex>(shift);
    }

    /// @returns link to all acceleration data structures - const access
    DETRAY_HOST_DEVICE constexpr auto accel_link() const
        -> const accel_link_type& {
        return m_accel_links;
    }

    /// @returns acc data structure link for a specific type of object - const
    template <ID obj_id>
    DETRAY_HOST_DEVICE constexpr auto accel_link() const -> const link_t& {
        return detail::get<obj_id>(m_accel_links);
    }

    /// Set surface finder link from @param link
    template <ID obj_id>
    DETRAY_HOST constexpr auto set_accel_link(const link_t link) -> void {
        detail::get<obj_id>(m_accel_links) = link;
    }

    /// Set surface finder link from @param id and @param index of the
    /// acceleration data structure (e.g. type and index of grid in surface
    /// store)
    template <ID obj_id>
    DETRAY_HOST constexpr auto set_accel_link(
        const typename link_t::id_type id,
        const typename link_t::index_type index) -> void {
        detail::get<obj_id>(m_accel_links) = link_t{id, index};
    }

    /// Equality operator
    ///
    /// @param rhs is the right-hand side to compare against.
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const volume_descriptor& rhs) const -> bool {
        return (m_id == rhs.m_id && m_index == rhs.m_index &&
                m_accel_links == rhs.m_accel_links);
    }

    private:
    /// How to interpret the boundary values
    volume_id m_id{volume_id::e_unknown};

    /// Volume index in the detector's volume container
    dindex m_index{dindex_invalid};

    /// Volume index in the detector's volume container
    dindex m_transform{dindex_invalid};

    /// Index range for every object type
    sf_link_type m_sf_links{};

    /// Links for every object type to an acceleration data structure
    accel_link_type m_accel_links{};

    /// Volume material
    detray::material<scalar_t> m_vol_mat = detray::vacuum<scalar_t>();
};

}  // namespace detray
