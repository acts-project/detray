
/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/barcode.hpp"
#include "detray/geometry/detail/surface_kernels.hpp"

namespace detray {

/// @brief Facade for a detray detector surface.
///
/// Provides an interface to geometry specific functionality like
/// local-to-global coordinate transforms or mask and material visitors. It
/// wraps a detector instance that contains the data and a surface descriptor
/// that contains the indices into the detector data containers for the
/// specific surface instance.
template <typename detector_t>  // @TODO: This needs a concept
class surface {

    /// Surface descriptor type
    using descr_t = typename detector_t::surface_type;

    using kernels = detail::surface_kernels<typename detector_t::transform3>;
    /// Vector type for track parameters in global coordinates
    using free_vector_type = typename kernels::free_vector_type;
    /// Vector type for track parameters in local (bound) coordinates
    using bound_vector_type = typename kernels::bound_vector_type;

    public:
    using algebra = typename detector_t::transform3;
    using transform3 = algebra;
    using point3 = typename transform3::point3;
    using vector3 = typename transform3::point3;
    using context = typename detector_t::geometry_context;

    /// Not allowed: always needs a detector and a descriptor.
    surface() = delete;

    /// Constructor from detector @param det and volume descriptor
    /// @param vol_idx from that detector.
    DETRAY_HOST_DEVICE
    constexpr surface(const detector_t &det, const descr_t &desc)
        : m_detector{det}, m_desc{desc} {}

    /// Constructor from detector @param det and volume index @param vol_idx in
    /// that detector.
    DETRAY_HOST_DEVICE
    constexpr surface(const detector_t &det, const geometry::barcode bcd)
        : surface(det, det.surface(bcd)) {}

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const surface &rhs) const -> bool {
        return (&m_detector == &(rhs.m_detector) and m_desc == rhs.m_desc);
    }

    /// @returns the surface barcode
    DETRAY_HOST_DEVICE
    constexpr auto barcode() const -> geometry::barcode {
        return m_desc.barcode();
    }

    /// @returns the index of the mother volume
    DETRAY_HOST_DEVICE
    constexpr auto volume() const -> dindex { return barcode().volume(); }

    /// @returns the index of the surface in the detector surface lookup
    DETRAY_HOST_DEVICE
    constexpr auto index() const -> dindex { return barcode().index(); }

    /// @returns the surface id (sensitive, passive or portal)
    DETRAY_HOST_DEVICE
    constexpr auto id() const -> surface_id { return barcode().id(); }

    /// @returns an id for the surface type (e.g. 'rectangle')
    DETRAY_HOST_DEVICE
    constexpr auto shape_id() const { return m_desc.mask().id(); }

    /// @returns the surface source link
    DETRAY_HOST_DEVICE
    constexpr auto source() const { return m_desc.source(); }

    /// @returns true if the surface is a senstive detector module.
    DETRAY_HOST_DEVICE
    constexpr auto is_sensitive() const -> bool {
        return barcode().id() == surface_id::e_sensitive;
    }

    /// @returns true if the surface is a portal.
    DETRAY_HOST_DEVICE
    constexpr auto is_portal() const -> bool {
        return barcode().id() == surface_id::e_portal;
    }

    /// @returns true if the surface is a passive detector element.
    DETRAY_HOST_DEVICE
    constexpr auto is_passive() const -> bool {
        return barcode().id() == surface_id::e_passive;
    }

    /// @returns the coordinate transform matrix of the surface
    DETRAY_HOST_DEVICE
    constexpr auto transform(const context &ctx) const -> const transform3 & {
        return m_detector.transform_store(ctx)[m_desc.transform()];
    }

    /// @returns the center position on the surface in global coordinates
    DETRAY_HOST_DEVICE
    constexpr auto center(const context &ctx) const -> const point3 & {
        return transform(ctx).translation();
    }

    /// @returns the surface normal in global coordinates
    DETRAY_HOST_DEVICE
    constexpr auto normal(const context &ctx) const -> const vector3 & {
        return transform(ctx).translation();
    }

    /// @returns the local position to the global point @param global for
    /// a given geometry context @param ctx and track direction @param dir
    DETRAY_HOST_DEVICE
    constexpr point3 global_to_local(const context &ctx, const point3 &global,
                                     const vector3 &dir) const {
        return visit_mask<typename kernels::global_to_local>(transform(ctx),
                                                             global, dir);
    }

    /// @returns the global position to the given local position @param local
    /// for a given geometry context @param ctx
    DETRAY_HOST_DEVICE
    constexpr point3 local_to_global(const context &ctx,
                                     const point3 &local) const {
        return visit_mask<typename kernels::local_to_global>(transform(ctx),
                                                             local);
    }

    /// @returns the track parametrization projected onto the surface (bound)
    DETRAY_HOST_DEVICE
    constexpr auto free_to_bound_vector(
        const context &ctx, const free_vector_type &free_vec) const {
        return visit_mask<typename kernels::free_to_bound_vector>(
            transform(ctx), free_vec);
    }

    /// @returns the global track parametrization from a bound representation
    DETRAY_HOST_DEVICE
    constexpr auto bound_to_free_vector(
        const context &ctx, const bound_vector_type &bound_vec) const {
        return visit_mask<typename kernels::bound_to_free_vector>(
            transform(ctx), bound_vec);
    }

    /// @returns the jacobian to go from a free to a bound track parametrization
    DETRAY_HOST_DEVICE
    constexpr auto free_to_bound_jacobian(
        const context &ctx, const free_vector_type &free_vec) const {
        return this
            ->template visit_mask<typename kernels::free_to_bound_jacobian>(
                transform(ctx), free_vec);
    }

    /// @returns the jacobian to go from a bound to a free track parametrization
    DETRAY_HOST_DEVICE
    constexpr auto bound_to_free_jacobian(
        const context &ctx, const bound_vector_type &bound_vec) const {
        return this
            ->template visit_mask<typename kernels::bound_to_free_jacobian>(
                transform(ctx), bound_vec);
    }

    /// @returns the path correction term
    DETRAY_HOST_DEVICE
    constexpr auto path_correction(const context &ctx, const vector3 &pos,
                                   const vector3 &dir,
                                   const vector3 &dtds) const {
        return visit_mask<typename kernels::path_correction>(transform(ctx),
                                                             pos, dir, dtds);
    }

    /// @returns the vertices in local frame
    DETRAY_HOST
    constexpr auto local_vertices() const {
        return visit_mask<typename kernels::local_vertices>();
    }

    /// @returns the vertices in global frame
    DETRAY_HOST
    constexpr auto global_vertices(const context &ctx) const {
        auto vertices = local_vertices();
        for (size_t i = 0; i < vertices.size(); i++){
            vertices[i] = local_to_global(ctx, vertices[i]);
        }
        return vertices;
    }

    /// Call a functor on the surfaces mask with additional arguments.
    ///
    /// @tparam functor_t the prescription to be applied to the mask
    /// @tparam Args      types of additional arguments to the functor
    template <typename functor_t, typename... Args>
    DETRAY_HOST_DEVICE constexpr auto visit_mask(Args &&... args) const {
        const auto &masks = m_detector.mask_store();

        return masks.template visit<functor_t>(m_desc.mask(),
                                               std::forward<Args>(args)...);
    }

    /// Call a functor on the surfaces material with additional arguments.
    ///
    /// @tparam functor_t the prescription to be applied to the mask
    /// @tparam Args      types of additional arguments to the functor
    template <typename functor_t, typename... Args>
    DETRAY_HOST_DEVICE constexpr auto visit_material(Args &&... args) const {
        const auto &materials = m_detector.material_store();

        return materials.template visit<functor_t>(m_desc.material(),
                                                   std::forward<Args>(args)...);
    }

    /// @returns a string stream that prints the surface details
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &os, const surface sf) {
        os << sf.barcode();
        os << " | trf.: " << sf.m_desc.transform();
        os << " | mask: " << static_cast<dindex>(sf.m_desc.mask().id()) << ", "
           << sf.m_desc.mask().index();
        os << " | mat.: " << static_cast<dindex>(sf.m_desc.material().id())
           << ", " << sf.m_desc.material().index();
        return os;
    }

    private:
    /// Access to the detector stores
    const detector_t &m_detector;
    /// Access to the descriptor
    const descr_t &m_desc;
};

template <typename detector_t, typename descr_t>
DETRAY_HOST_DEVICE surface(const detector_t &, const descr_t &)
    ->surface<detector_t>;

template <typename detector_t>
DETRAY_HOST_DEVICE surface(const detector_t &, const geometry::barcode)
    ->surface<detector_t>;

}  // namespace detray
