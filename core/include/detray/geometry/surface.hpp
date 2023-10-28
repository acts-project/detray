
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

// System include(s)
#include <ostream>
#include <type_traits>

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
    using scalar_type = typename detector_t::scalar_type;
    using transform3 = algebra;
    using point2 = typename algebra::point2;
    using point3 = typename algebra::point3;
    using vector3 = typename algebra::point3;
    using context = typename detector_t::geometry_context;

    /// Not allowed: always needs a detector and a descriptor.
    surface() = delete;

    /// Constructor from detector @param det and surface descriptor
    /// @param desc from that detector.
    DETRAY_HOST_DEVICE
    constexpr surface(const detector_t &det, const descr_t &desc)
        : m_detector{det}, m_desc{desc} {}

    /// Constructor from detector @param det and barcode @param bcd in
    /// that detector.
    DETRAY_HOST_DEVICE
    constexpr surface(const detector_t &det, const geometry::barcode bcd)
        : surface(det, det.surface(bcd)) {}

    /// Conversion to surface interface around constant detector type
    template <typename detector_type = detector_t,
              std::enable_if_t<!std::is_const_v<detector_type>, bool> = true>
    DETRAY_HOST_DEVICE constexpr operator surface<const detector_type>() const {
        return surface<const detector_type>{this->m_detector, this->m_desc};
    }

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

    /// @returns the extra bits in the barcode
    DETRAY_HOST_DEVICE
    constexpr auto extra() const -> dindex { return barcode().extra(); }

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

    /// @returns the mask volume link
    DETRAY_HOST_DEVICE
    constexpr auto volume_link() const {
        return visit_mask<typename kernels::get_volume_link>();
    }

    /// @returns the coordinate transform matrix of the surface
    DETRAY_HOST_DEVICE
    constexpr auto transform(const context &ctx) const -> const transform3 & {
        return m_detector.transform_store(ctx)[m_desc.transform()];
    }

    /// @returns the center position on the surface in global coordinates
    DETRAY_HOST_DEVICE
    constexpr auto center(const context &ctx) const -> const point3 {
        return transform(ctx).translation();
    }

    /// @returns the surface normal in global coordinates at a given bound/local
    /// position @param p
    template <typename point_t = point2,
              std::enable_if_t<std::is_same_v<point_t, point3> or
                                   std::is_same_v<point_t, point2>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr auto normal(const context &ctx,
                                             const point_t &p) const
        -> const vector3 {
        return visit_mask<typename kernels::normal>(transform(ctx), p);
    }

    /// @returns the cosine of the incidence angle given a local/bound position
    /// @param p and a global direction @param dir
    template <typename point_t = point2,
              std::enable_if_t<std::is_same_v<point_t, point3> or
                                   std::is_same_v<point_t, point2>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr auto cos_angle(const context &ctx,
                                                const vector3 &dir,
                                                const point_t &p) const
        -> scalar_type {
        return vector::dot(dir, normal(ctx, p));
    }

    /// @returns the bound (2D) position to the global point @param global for
    /// a given geometry context @param ctx and track direction @param dir
    DETRAY_HOST_DEVICE
    constexpr point2 global_to_bound(const context &ctx, const point3 &global,
                                     const vector3 &dir) const {
        return visit_mask<typename kernels::global_to_bound>(transform(ctx),
                                                             global, dir);
    }

    /// @returns the local position to the global point @param global for
    /// a given geometry context @param ctx and track direction @param dir
    DETRAY_HOST_DEVICE
    constexpr point3 global_to_local(const context &ctx, const point3 &global,
                                     const vector3 &dir) const {
        return visit_mask<typename kernels::global_to_local>(transform(ctx),
                                                             global, dir);
    }

    /// @returns the global position to the given local/bound position @param p
    /// for a given geometry context @param ctx
    template <typename point_t,
              std::enable_if_t<std::is_same_v<point_t, point3> or
                                   std::is_same_v<point_t, point2>,
                               bool> = true>
    DETRAY_HOST_DEVICE constexpr point3 local_to_global(
        const context &ctx, const point_t &p, const vector3 &dir) const {
        return visit_mask<typename kernels::local_to_global>(transform(ctx), p,
                                                             dir);
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

    /// @returns the vertices in local frame with @param n_seg the number of
    /// segments used along acrs
    DETRAY_HOST
    constexpr auto local_vertices(const dindex n_seg) const {
        return visit_mask<typename kernels::vertices>(n_seg);
    }

    /// @returns the vertices in global frame with @param n_seg the number of
    /// segments used along acrs
    DETRAY_HOST
    constexpr auto global_vertices(const context &ctx, const dindex n_seg,
                                   const vector3 &dir) const {
        auto vertices = local_vertices(n_seg);
        for (size_t i = 0; i < vertices.size(); i++) {
            vertices[i] = local_to_global(ctx, vertices[i], dir);
        }
        return vertices;
    }

    /// @brief Lower and upper point for minimal axis aligned bounding box.
    ///
    /// Computes the min and max vertices in a local cartesian frame.
    DETRAY_HOST
    constexpr auto local_min_bounds(
        const scalar_type env =
            std::numeric_limits<scalar_type>::epsilon()) const {
        return visit_mask<typename kernels::local_min_bounds>(env);
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

    /// Do a consistency check on the surface after building the detector.
    ///
    /// @param os output stream for error messages.
    ///
    /// @returns true if the surface is consistent
    DETRAY_HOST bool self_check(std::ostream &os) const {
        if (barcode().is_invalid()) {
            os << "ERROR: Invalid barcode for surface:\n" << *this << std::endl;
            return false;
        }
        if (index() >= m_detector.surface_lookup().size()) {
            os << "ERROR: Surface index out of bounds for surface:\n"
               << *this << std::endl;
            return false;
        }
        if (volume() >= m_detector.volumes().size()) {
            os << "ERROR: Surface volume index out of bounds for surface:\n"
               << *this << std::endl;
            return false;
        }
        if (is_invalid_value(m_desc.transform())) {
            os << "ERROR: Surface transform undefined for surface:\n"
               << *this << std::endl;
            return false;
        }
        if (m_desc.transform() >= m_detector.transform_store().size()) {
            os << "ERROR: Surface transform index out of bounds for surface:\n"
               << *this << std::endl;
            return false;
        }
        if (is_invalid_value(m_desc.mask())) {
            os << "ERROR: Surface does not have a valid mask link:\n"
               << *this << std::endl;
            return false;
        }
        // Only check, if there is material in the detector
        if (not m_detector.material_store().all_empty()) {
            if (is_invalid_value(m_desc.material())) {
                os << "ERROR: Surface does not have valid material:\n"
                   << *this << std::endl;
                return false;
            }
        }
        // Check the mask boundaries
        if (not visit_mask<typename kernels::mask_self_check>(os)) {
            os << "\nSurface: " << *this << std::endl;
            return false;
        }
        // Check the mask volume link
        const auto vol_link = visit_mask<typename kernels::get_volume_link>();
        if (is_portal()) {
            if (vol_link == volume()) {
                os << "ERROR: Portal surface links to mother volume:\n"
                   << *this << std::endl;
                return false;
            }
        } else if (vol_link != volume()) {
            os << "ERROR: Passive/sensitive surface does not link to mother "
                  "volume:\n"
               << *this << std::endl;
            return false;
        }

        return true;
    }

    /// @returns a string stream that prints the surface details
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &os, const surface &sf) {
        os << sf.m_desc;
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
