/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/detail/volume_kernels.hpp"
#include "detray/materials/material.hpp"

// System include(s)
#include <iostream>
#include <sstream>

namespace detray {

/// @brief Facade for a detray detector volume.
///
/// Volume class that acts as a logical container in the detector for geometry
/// objects, i.e. surfaces. The volume boundary surfaces, the so-called portals,
/// carry index links that join adjacent volumes. The volume class itself does
/// not contain any data itself, but keeps a descriptor with index-based links
/// into the data containers that are managed by the detector.
/// Every type of surface that is known by the volume (determined by the type
/// and size of the ID enum in the descriptor) lives in it's own geometry
/// accelerator data structure, e.g. portals reside in a brute force
/// accelerator (a simple vector), while sensitive surfaces are usually sorted
/// into a spacial grid.
///
/// @TODO: Add access to the volume placement transform and volume center
template <typename detector_t>  // @TODO: This needs a concept
class detector_volume {

    /// Scalar type
    using scalar_type = typename detector_t::scalar_type;

    /// Volume descriptor type
    using descr_t = typename detector_t::volume_type;

    public:
    /// In case the geometry needs to be printed
    using name_map = dmap<dindex, std::string>;

    /// Allow detector to access descriptor. @TODO: Remove once possible
    friend detector_t;

    /// Not allowed: always needs a detector and a descriptor.
    detector_volume() = delete;

    /// Constructor from detector @param det and volume descriptor
    /// @param vol_idx from that detector.
    constexpr detector_volume(const detector_t &det, const descr_t &desc)
        : m_detector{det}, m_desc{desc} {}

    /// Constructor from detector @param det and volume index @param vol_idx in
    /// that detector.
    constexpr detector_volume(const detector_t &det, const dindex vol_idx)
        : detector_volume(det, det.volumes()[vol_idx]) {}

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE
    constexpr auto operator==(const detector_volume &rhs) const -> bool {
        return (&m_detector == &(rhs.m_detector) and m_desc == rhs.m_desc);
    }

    /// @returns the volume shape id, e.g. 'cylinder'.
    DETRAY_HOST_DEVICE
    constexpr auto id() const -> volume_id { return m_desc.id(); }

    /// @returns the index of the volume in the detector volume container.
    DETRAY_HOST_DEVICE
    constexpr auto index() const -> dindex { return m_desc.index(); }

    /// @returns the volume name (add an offset for the detector name).
    DETRAY_HOST_DEVICE
    constexpr auto name(const name_map &names) const -> const std::string & {
        return names.at(m_desc.index() + 1u);
    }

    /// @returns the (non contextual) transform for the placement of the
    /// volume in the detector geometry.
    DETRAY_HOST_DEVICE
    constexpr auto transform() const -> const
        typename detector_t::transform3 & {
        return m_detector.transform_store({})[m_desc.transform()];
    }

    /// @returns the center point of the volume.
    DETRAY_HOST_DEVICE
    constexpr auto center() const -> typename detector_t::point3 {
        return transform().translation();
    }

    /// @returns the material of the volume
    DETRAY_HOST_DEVICE
    constexpr auto material() const -> material<scalar_type> {
        return m_desc.material();
    }

    /// Apply a functor to all surfaces in the volume.
    ///
    /// @tparam functor_t the prescription to be applied to the surfaces
    /// @tparam Args      types of additional arguments to the functor
    template <typename functor_t,
              int I = static_cast<int>(descr_t::object_id::e_size) - 1,
              typename... Args>
    DETRAY_HOST_DEVICE constexpr void visit_surfaces(Args &&... args) const {
        visit_surfaces_impl<detail::surface_getter<functor_t>>(
            std::forward<Args>(args)...);
    }

    /// Apply a functor to a neighborhood of surfaces around a track position
    /// in the volume.
    ///
    /// @tparam functor_t the prescription to be applied to the surfaces (
    ///                   customization point for the navigation)
    /// @tparam track_t   the track around which to build up the neighborhood
    /// @tparam Args      types of additional arguments to the functor
    template <typename functor_t,
              int I = static_cast<int>(descr_t::object_id::e_size) - 1,
              typename track_t, typename... Args>
    DETRAY_HOST_DEVICE constexpr void visit_neighborhood(
        const track_t &track, Args &&... args) const {
        visit_surfaces_impl<detail::neighborhood_getter<functor_t>>(
            m_detector, m_desc, track, std::forward<Args>(args)...);
    }

    /// @returns the maximum number of surface candidates during a neighborhood
    /// lookup
    // TODO: Remove
    template <int I = static_cast<int>(descr_t::object_id::e_size) - 1>
    DETRAY_HOST_DEVICE constexpr auto n_max_candidates(
        unsigned int n = 0u) const -> unsigned int {
        // Get the index of the surface collection with type index 'I'
        constexpr auto sf_col_id{
            static_cast<typename detector_t::accel::id>(I)};
        const auto &link{
            m_desc
                .template link<static_cast<typename descr_t::object_id>(I)>()};

        // Check if this volume holds such a collection and, if so, add max
        // number of candidates that we can expect from it
        if (not link.is_invalid()) {
            const unsigned int n_max{
                m_detector.accelerator_store()
                    .template get<sf_col_id>()[detail::get<1>(link)]
                    .n_max_candidates()};
            // @todo: Remove when local navigation becomes available !!!!
            n += n_max > 20u ? 20u : n_max;
        }
        // Check the next surface collection type
        if constexpr (I > 0) {
            return n_max_candidates<I - 1>(n);
        } else {
            return n;
        }
    }

    /// Do a consistency check on the volume after building the detector.
    ///
    /// @param os output stream for error messages.
    ///
    /// @returns true if the volume is consistent
    DETRAY_HOST bool self_check(std::ostream &os) const {
        if (id() == volume_id::e_unknown) {
            os << "ERROR: Unknown volume shape type in volume:\n"
               << *this << std::endl;
            return false;
        }
        if (detail::is_invalid_value(index())) {
            os << "ERROR: Volume index undefined in volume:\n"
               << *this << std::endl;
            return false;
        }
        if (index() >= m_detector.volumes().size()) {
            os << "ERROR: Volume index out of bounds in volume:\n"
               << *this << std::endl;
            return false;
        }
        if (detail::is_invalid_value(m_desc.transform())) {
            os << "ERROR: Volume transform undefined in volume:\n"
               << *this << std::endl;
            return false;
        }
        if (m_desc.transform() >= m_detector.transform_store().size()) {
            os << "ERROR: Volume transform index out of bounds in volume:\n"
               << *this << std::endl;
            return false;
        }
        const auto &acc_link = m_desc.full_link();
        if (detail::is_invalid_value(acc_link[0])) {
            os << "ERROR: Link to portal lookup broken in volume: " << acc_link
               << "\n in volume: " << *this << std::endl;
            return false;
        }

        // Warnings
        bool suspicious_links = false;
        std::stringstream warnigns{};
        for (std::size_t i = 1u; i < acc_link.size(); ++i) {
            // An acceleration data structure link was set, but is invalid
            if (!acc_link[i].is_invalid_id() and
                acc_link[i].is_invalid_index()) {
                suspicious_links = true;
                warnigns << "Link to acceleration data structure "
                         << static_cast<int>(acc_link[i].id()) << " is invalid"
                         << std::endl;
            }
        }
        if (suspicious_links) {
            std::cout << "WARNING: " << warnigns.str()
                      << " in volume: " << *this << std::endl;
        }

        return true;
    }

    /// @returns a string stream that prints the volume details
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &os,
                                    const detector_volume &v) {
        os << "id: " << static_cast<int>(v.m_desc.id());
        os << " | index: " << v.m_desc.index();
        os << " | trf.: " << v.m_desc.transform();
        os << " | acc link: " << v.m_desc.full_link();

        return os;
    }

    private:
    /// Apply a functor to all acceleration structures of this volume.
    ///
    /// @tparam functor_t the prescription to be applied to the acc structure
    /// @tparam Args      types of additional arguments to the functor
    template <typename functor_t,
              int I = static_cast<int>(descr_t::object_id::e_size) - 1,
              typename... Args>
    DETRAY_HOST_DEVICE constexpr void visit_surfaces_impl(
        Args &&... args) const {
        // Get the acceleration data structures for this volume
        const auto &link{
            m_desc
                .template link<static_cast<typename descr_t::object_id>(I)>()};

        // Only visit, if object type is contained in volume
        if (not link.is_invalid()) {
            // Run over the surfaces in a single acceleration data structure
            // and apply the functor to the resulting neighborhood
            m_detector.accelerator_store().template visit<functor_t>(
                link, std::forward<Args>(args)...);
        }
        // Check the next surface type
        if constexpr (I > 0) {
            visit_surfaces_impl<functor_t, I - 1, Args...>(
                std::forward<Args>(args)...);
        }
    }

    /// Access to the detector stores
    const detector_t &m_detector;
    /// Access to the descriptor
    const descr_t &m_desc;
};

}  // namespace detray
