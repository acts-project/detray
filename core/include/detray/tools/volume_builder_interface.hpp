/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/geometry/detector_volume.hpp"

// System include(s)
#include <memory>

namespace detray {

template <typename detector_t>
class surface_factory_interface;

/// @brief Interface for volume builders (and volume builder decorators)
template <typename detector_t>
class volume_builder_interface {

    public:
    using scalar_type = typename detector_t::scalar_type;
    template <typename T, std::size_t N>
    using array_type = typename detector_t::template array_type<T, N>;

    virtual ~volume_builder_interface() = default;

    /// @brief Initializes a new volume with shape id @param id in detector
    /// @param det
    DETRAY_HOST
    virtual void init_vol(detector_t &det, const volume_id id) = 0;

    /// @returns the global index for the volume
    /// @note the correct index is only available after calling @c init_vol
    DETRAY_HOST
    virtual auto get_vol_index() -> dindex = 0;

    /// @returns reading access to the volume
    DETRAY_HOST
    virtual auto operator()() const -> const
        typename detector_t::volume_type & = 0;
    DETRAY_HOST
    virtual auto operator()() -> typename detector_t::volume_type & = 0;

    /// @brief Adds a volume and all of its contents to a detector
    DETRAY_HOST
    virtual auto build(detector_t &det,
                       typename detector_t::geometry_context ctx = {}) ->
        typename detector_t::volume_type * = 0;

    /// @brief Add the transform for the volume placement - copy
    DETRAY_HOST
    virtual void add_volume_placement(
        const typename detector_t::transform3 &trf = {}) = 0;

    /// @brief Add the transform for the volume placement from the translation
    /// @param t. The rotation will be the identity matrix.
    DETRAY_HOST
    virtual void add_volume_placement(const typename detector_t::point3 &t) = 0;

    /// @brief Add the transform for the volume placement from the translation
    /// @param t , the new x- and z-axis (@param x, @param z).
    DETRAY_HOST
    virtual void add_volume_placement(
        const typename detector_t::point3 &t,
        const typename detector_t::vector3 &x,
        const typename detector_t::vector3 &z) = 0;

    /// @brief Add the portals to the volume.
    /// @returns the index range of the portals in the temporary surface
    /// container used by the factory ( gets final update in @c build() )
    DETRAY_HOST
    virtual void add_portals(
        std::shared_ptr<surface_factory_interface<detector_t>> pt_factory,
        typename detector_t::geometry_context ctx = {}) = 0;

    /// @brief Add sensitive surfaces to the volume, if any
    /// @returns the index range of the sensitives in the temporary surface
    /// container used by the factory ( gets final update in @c build() )
    DETRAY_HOST
    virtual void add_sensitives(
        std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
        typename detector_t::geometry_context ctx = {}) = 0;

    /// @brief Add passive surfaces to the volume, if any
    /// @returns the index range of the passives in the temporary surface
    /// container used by the factory ( gets final update in @c build() )
    DETRAY_HOST
    virtual void add_passives(
        std::shared_ptr<surface_factory_interface<detector_t>> ps_factory,
        typename detector_t::geometry_context ctx = {}) = 0;
};

/// @brief Decorator for the volume builder.
///
/// Can be volume builders that introduce special sorting/memory layout, or
/// accelerator builders, like the grid builder.
template <typename detector_t>
class volume_decorator : public volume_builder_interface<detector_t> {

    public:
    using scalar_type = typename detector_t::scalar_type;
    template <typename T, std::size_t N>
    using array_type = typename detector_t::template array_type<T, N>;

    DETRAY_HOST
    volume_decorator(
        std::unique_ptr<volume_builder_interface<detector_t>> vol_builder)
        : m_builder(std::move(vol_builder)) {}

    /// Overwrite interface functions using callbacks to the volume builder
    /// @{
    DETRAY_HOST
    void init_vol(detector_t &det, const volume_id id) override {
        m_builder->init_vol(det, id);
    }

    DETRAY_HOST
    auto operator()() -> typename detector_t::volume_type & override {
        return m_builder->operator()();
    }

    DETRAY_HOST
    auto operator()() const -> const
        typename detector_t::volume_type & override {
        return m_builder->operator()();
    }

    DETRAY_HOST
    auto get_vol_index() -> dindex override {
        return m_builder->get_vol_index();
    }

    DETRAY_HOST
    auto build(detector_t &det,
               typename detector_t::geometry_context /*ctx*/ = {}) ->
        typename detector_t::volume_type * override {
        return m_builder->build(det);
    }

    DETRAY_HOST
    void add_volume_placement(
        const typename detector_t::transform3 &trf = {}) override {
        return m_builder->add_volume_placement(trf);
    }

    DETRAY_HOST
    void add_volume_placement(const typename detector_t::point3 &t) override {
        return m_builder->add_volume_placement(t);
    }

    DETRAY_HOST
    void add_volume_placement(const typename detector_t::point3 &t,
                              const typename detector_t::vector3 &x,
                              const typename detector_t::vector3 &z) override {
        return m_builder->add_volume_placement(t, x, z);
    }

    DETRAY_HOST
    void add_portals(
        std::shared_ptr<surface_factory_interface<detector_t>> pt_factory,
        typename detector_t::geometry_context ctx = {}) override {
        return m_builder->add_portals(std::move(pt_factory), ctx);
    }

    DETRAY_HOST
    void add_sensitives(
        std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
        typename detector_t::geometry_context ctx = {}) override {
        return m_builder->add_sensitives(std::move(sf_factory), ctx);
    }

    DETRAY_HOST
    void add_passives(
        std::shared_ptr<surface_factory_interface<detector_t>> ps_factory,
        typename detector_t::geometry_context ctx = {}) override {
        return m_builder->add_passives(std::move(ps_factory), ctx);
    }
    /// @}

    protected:
    std::unique_ptr<volume_builder_interface<detector_t>> m_builder;
};

namespace detail {

/// A functor to update the material index in surface objects
struct material_index_update {

    template <typename group_t, typename index_t, typename surface_t>
    DETRAY_HOST inline void operator()(const group_t &group,
                                       const index_t & /*index*/,
                                       surface_t &sf) const {
        sf.update_material(static_cast<dindex>(group.size()));
    }
};

}  // namespace detail

}  // namespace detray
