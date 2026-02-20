/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/navigation.hpp"
#include "detray/geometry/barcode.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/logging.hpp"

namespace detray::navigation {

/// Result type of a navigation update
template <typename metadata_t, std::size_t N = 1u>
struct alignas(128) result {
    using detector_t = detector<metadata_t, device_container_types>;
    using nav_link_t = typename detector_t::surface_type::navigation_link;
    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;

    /// TODO: replace with @c dsimd<algebra_t, T>
    template <typename T>
    using simd_t = std::array<T, N>;

    /// @returns current navigation status - const
    DETRAY_HOST_DEVICE
    constexpr auto status(const std::size_t i = 0u) const
        -> navigation::status {
        return m_status[i];
    }

    /// @returns the navigation heartbeat (is this stream still being updated?)
    DETRAY_HOST_DEVICE
    constexpr bool is_alive(const std::size_t i = 0u) const {
        return status(i) > navigation::status::e_unknown;
    }

    /// @returns the externally provided mask tolerance - const
    DETRAY_HOST_DEVICE
    constexpr scalar_t external_tol(const std::size_t i = 0u) const {
        return m_external_mask_tol[i];
    }

    /// Set externally provided mask tolerance according to noise prediction
    DETRAY_HOST_DEVICE
    constexpr void set_external_tol(const scalar_t tol,
                                    const std::size_t i = 0u) {
        DETRAY_VERBOSE_HOST("Setting external mask tolerance: " << tol);
        m_external_mask_tol[i] = tol;
    }
    /// @returns navigation trust level - const
    DETRAY_HOST_DEVICE
    constexpr auto trust_level(const std::size_t i = 0u) const
        -> navigation::trust_level {
        return m_trust_level[i];
    }

    /// Update navigation trust level to 'no trust'
    /// @{
    DETRAY_HOST_DEVICE
    constexpr void set_no_trust(const std::size_t i = 0u) {
        DETRAY_VERBOSE_HOST_DEVICE("Flagging re-navigation: \"no trust\"");
        m_trust_level[i] = navigation::trust_level::e_no_trust;
    }

    /// @note only two trust levels for basic navigation implementation
    DETRAY_HOST_DEVICE
    constexpr void set_high_trust(const std::size_t i = 0u) {
        m_trust_level[i] = navigation::trust_level::e_high;
    }

    /// @note only two trust levels for basic navigation implementation
    DETRAY_HOST_DEVICE
    constexpr void set_fair_trust(const std::size_t i = 0u) {
        m_trust_level[i] = navigation::trust_level::e_fair;
    }
    /// @}

    /// Helper method to check the track has reached a module surface
    DETRAY_HOST_DEVICE
    constexpr auto is_on_surface(const std::size_t i = 0u) const -> bool {
        return (m_status[i] == navigation::status::e_on_object ||
                m_status[i] == navigation::status::e_on_portal);
    }

    /// Helper method to check the track has reached a sensitive surface
    DETRAY_HOST_DEVICE
    constexpr auto is_on_sensitive(const std::size_t i = 0u) const -> bool {
        return (m_status[i] == navigation::status::e_on_object) &&
               (barcode(i).id() == surface_id::e_sensitive);
    }

    /// Helper method to check the track has reached a passive surface
    DETRAY_HOST_DEVICE
    constexpr auto is_on_passive(const std::size_t i = 0u) const -> bool {
        return (m_status[i] == navigation::status::e_on_object) &&
               (barcode(i).id() == surface_id::e_passive);
    }

    /// Helper method to check the track has reached a portal surface
    DETRAY_HOST_DEVICE
    constexpr auto is_on_portal(const std::size_t i = 0u) const -> bool {
        return m_status[i] == navigation::status::e_on_portal;
    }

    /// Helper method to check the track has encountered material
    DETRAY_HOST_DEVICE
    constexpr auto encountered_sf_material(const std::size_t i = 0u) const
        -> bool {
        return is_on_surface(i) && m_current[i].has_material();
    }

    /// @returns flag that indicates whether navigation was successful
    DETRAY_HOST_DEVICE
    constexpr bool finished(const std::size_t i = 0u) const {
        return (m_status[i] == navigation::status::e_exit);
    }

    /// @returns current volume index - const
    DETRAY_HOST_DEVICE
    constexpr auto volume(const std::size_t i = 0u) const -> nav_link_t {
        return m_volume_index[i];
    }

    /// @returns current/previous object that was reached - const
    DETRAY_HOST_DEVICE
    constexpr auto current(const std::size_t i = 0u) const {
        assert(is_on_surface());
        return m_current[i];
    }

    /// @returns next object that we want to reach (current target) - const
    DETRAY_HOST_DEVICE
    constexpr auto target(const std::size_t i = 0u) const -> geometry::barcode {
        return m_target[i];
    }

    /// @returns barcode of the detector surface the navigator is on
    /// (invalid when not on surface) - const
    DETRAY_HOST_DEVICE
    constexpr auto barcode(const std::size_t i = 0u) const
        -> geometry::barcode {
        return m_current[i].barcode();
    }

    /// @returns true if the current candidate lies on the surface edge
    DETRAY_HOST_DEVICE
    constexpr bool is_edge_candidate(const std::size_t i = 0u) const {
        assert(is_on_surface(i));
        return m_current_status[i] == intersection::status::e_edge;
    }

    /// @returns true if the current candidate lies on the surface
    DETRAY_HOST_DEVICE
    constexpr bool is_good_candidate(const std::size_t i = 0u) const {
        assert(is_on_surface());
        return m_current_status[i] == intersection::status::e_inside;
    }

    /// @returns true if the current candidate lies on the surface,
    /// inlcuding its edge
    DETRAY_HOST_DEVICE
    constexpr bool is_probably_candidate(const std::size_t i = 0u) const {
        assert(is_on_surface());
        return m_current_status[i] <= intersection::status::e_edge;
    }

    /// Scalar representation of the navigation state,
    /// @returns distance to next target
    DETRAY_HOST_DEVICE
    constexpr scalar_t operator()(const std::size_t i = 0u) const {
        assert(math::isfinite(m_dist_to_next[i]));
        return m_dist_to_next[i];
    }

    /// @param surface_t the surface interface type that is required
    /// @returns current detector surface the navigator is on - const
    template <template <typename> class surface_t = tracking_surface>
    DETRAY_HOST_DEVICE constexpr auto current_surface(
        const std::size_t i = 0u) const {
        assert(is_on_surface(i));
        return surface_t<detector_t>{*m_detector, m_current[i]};
    }

    /// @param volume_t the volume interface type that is required
    /// @returns current detector volume of the navigation stream
    template <template <typename> class volume_t = tracking_volume>
    DETRAY_HOST_DEVICE constexpr auto current_volume(
        const std::size_t i = 0u) const {
        return volume_t<detector_t>{*m_detector, m_volume_index[i]};
    }

    /// @param surface_t the surface interface type that is required
    /// @returns the next surface the navigator intends to reach
    template <template <typename> class surface_t = tracking_surface>
    DETRAY_HOST_DEVICE constexpr auto next_surface(
        const std::size_t i = 0u) const {
        return surface_t<detector_t>{*m_detector, m_target[i]};
    }

    /// Navigation reaches final target or leaves detector world. Stop
    /// navigation.
    DETRAY_HOST_DEVICE constexpr void exit(const std::size_t i = 0u) {
        m_status[i] = navigation::status::e_exit;
        DETRAY_VERBOSE_HOST_DEVICE("Exited");
    }

    /// Navigation is being paused by actor: Maintain the navigation state
    /// and resume later
    DETRAY_HOST_DEVICE constexpr void pause(const std::size_t i = 0u) {
        DETRAY_VERBOSE_HOST_DEVICE("Paused by actor");
    }

    /// Navigation state that cannot be recovered from. Leave the other
    /// data for inspection.
    ///
    /// @param custom_msg additional information on the reason for the error
    DETRAY_HOST_DEVICE constexpr void abort(
        const char *custom_msg = "Navigator (unkown reason)",
        const std::size_t i = 0u) {
        m_status[i] = navigation::status::e_abort;

        /// Wrapper around the custom message that a print inspector can
        /// understand
        struct message_wrapper {
            const char *const m_msg{nullptr};

            DETRAY_HOST_DEVICE
            constexpr const char *operator()() const { return m_msg; }
        };

        assert(custom_msg != nullptr);
        DETRAY_ERROR_HOST("Aborted: " << custom_msg);
    }

    /// Navigation state that cannot be recovered from. Leave the other
    /// data for inspection.
    ///
    /// @param debug_msg_generator functor that returns additional
    ///                            information on the reason for the error
    template <typename debug_msg_generator_t>
        requires(!std::same_as<char *, debug_msg_generator_t>)
    DETRAY_HOST_DEVICE constexpr void abort(
        const debug_msg_generator_t &debug_msg_generator,
        const std::size_t i = 0u) {
        m_status[i] = navigation::status::e_abort;

        DETRAY_ERROR_HOST("Aborted: " << debug_msg_generator());
    }

    simd_t<typename detector_t::surface_type> m_current{};
    simd_t<geometry::barcode> m_target{};
    simd_t<scalar_t> m_dist_to_next{detail::invalid_value<scalar_t>()};
    simd_t<scalar_t> m_external_mask_tol{0.f * unit<scalar_t>::mm};
    simd_t<nav_link_t> m_volume_index{detail::invalid_value<nav_link_t>()};
    simd_t<intersection::status> m_current_status{
        intersection::status::e_outside};
    simd_t<navigation::status> m_status{navigation::status::e_unknown};
    simd_t<navigation::trust_level> m_trust_level{
        navigation::trust_level::e_no_trust};
    detector_t *m_detector{nullptr};
};

}  // namespace detray::navigation
