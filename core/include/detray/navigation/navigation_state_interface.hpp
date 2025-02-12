/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/barcode.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/utils/ranges.hpp"

namespace detray {

/// @brief Facade for a navigator state towards the propagation flow and actors.
///
/// The state interface extends the base navigation state implementation, which
/// contains the state data of a single navigation stream. It is accessible to
/// the actors in the propagation, for which it defines the public interface
/// towards the navigation state.
/// The navigation state interface also provides an iterable view on the
/// reachable candidates in the navigation state (only const access).
///
/// @tparam navigation_state_t data state base type of a navigator state
template <class navigation_state_t>
class navigation_state_interface
    : public detray::ranges::view_interface<
          navigation_state_interface<navigation_state_t>> {

    /// Cast to the navigation state type to access its methods
    /// @{
    DETRAY_HOST_DEVICE
    constexpr auto cast_impl() -> navigation_state_t& {
        return static_cast<navigation_state_t&>(*this);
    }
    DETRAY_HOST_DEVICE
    constexpr auto cast_impl() const -> const navigation_state_t& {
        return static_cast<const navigation_state_t&>(*this);
    }
    /// @}

    public:
    /// @returns a reference to the detector that is being navigated - const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) detector() const {
        return cast_impl().get_detector();
    }

    /// @returns the navigation heartbeat (is this stream still being updated?)
    DETRAY_HOST_DEVICE
    constexpr bool is_alive() const { return cast_impl().get_heartbeat(); }

    /// @returns the navigation status (on object, towards object etc)
    DETRAY_HOST_DEVICE
    constexpr auto status() const { return cast_impl().get_status(); }

    /// @returns the navigation trust level (how outdated is the state?)
    DETRAY_HOST_DEVICE
    constexpr auto trust_level() const { return cast_impl().get_trust_level(); }

    /// @returns current navigation direction - const
    DETRAY_HOST_DEVICE
    constexpr auto direction() const -> navigation::direction {
        return cast_impl().get_direction();
    }

    /// Change direction in which the navigator should search for candidates
    DETRAY_HOST_DEVICE
    constexpr void direction(const navigation::direction dir) {
        cast_impl().set_direction(dir);
    }

    /// @returns currently cached candidates - const
    DETRAY_HOST_DEVICE
    constexpr const auto& candidates() const {
        return cast_impl().get_candidates();
    }

    /// @returns number of candidates the cache can hold - const
    DETRAY_HOST_DEVICE
    constexpr dindex n_candidates() const {
        return cast_impl().get_n_candidates();
    }

    /// @returns true if there are no cached candidates left - const
    DETRAY_HOST_DEVICE
    constexpr bool is_exhausted() const { return n_candidates() == 0u; }

    /// @returns start position of the valid candidate range - const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) begin() const { return cast_impl().get_begin(); }

    /// @returns sentinel of the valid candidate range - const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) end() const { return cast_impl().get_end(); }

    /// @returns current volume index - const
    DETRAY_HOST_DEVICE
    constexpr auto volume() const { return cast_impl().get_volume_idx(); }

    /// Set start/new volume from volume index @param v
    DETRAY_HOST_DEVICE
    constexpr void set_volume(dindex v) { cast_impl().set_volume_idx(v); }

    /// @returns current candidate (intersection) that was reached - const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) current() const {
        assert(is_on_surface());
        return cast_impl().get_current();
    }

    /// @returns next candidate (intersection) that we want to reach - const
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) target() const { return cast_impl().get_target(); }

    /// Scalar representation of the navigation state,
    /// @returns distance to target
    DETRAY_HOST_DEVICE
    constexpr auto operator()() const {
        return static_cast<float>(direction()) * target().path;
    }

    /// @returns barcode of the detector surface the navigator is on
    /// (invalid when not on surface)
    DETRAY_HOST_DEVICE
    constexpr auto barcode() const -> geometry::barcode {
        return current().sf_desc.barcode();
    }

    /// Lower the navigation trust level to 'no trust'
    DETRAY_HOST_DEVICE
    constexpr void set_no_trust() {
        cast_impl().set_trust_level(navigation::trust_level::e_no_trust);
    }

    /// Lower the navigation trust level to 'high trust'
    DETRAY_HOST_DEVICE
    constexpr void set_high_trust() {
        cast_impl().set_trust_level(trust_level() <
                                            navigation::trust_level::e_high
                                        ? trust_level()
                                        : navigation::trust_level::e_high);
    }

    /// Lower the navigation trust level to 'fair trust'
    DETRAY_HOST_DEVICE
    constexpr void set_fair_trust() {
        cast_impl().set_trust_level(trust_level() <
                                            navigation::trust_level::e_fair
                                        ? trust_level()
                                        : navigation::trust_level::e_fair);
    }

    /// @returns @c true if the track has reached a detector surface
    DETRAY_HOST_DEVICE
    constexpr bool is_on_surface() const { return cast_impl().on_surface(); }

    /// @returns @c true if the track has reached a sensitive surface
    DETRAY_HOST_DEVICE
    constexpr bool is_on_sensitive() const {
        return cast_impl().on_surface() &&
               cast_impl().get_current().sf_desc.is_sensitive();
    }

    /// @returns @c true if the track has reached a passive surface
    DETRAY_HOST_DEVICE
    constexpr bool is_on_passive() const {
        return cast_impl().on_surface() && current().sf_desc.is_passive();
    }

    /// @returns @c true if the track has reached a portal surface
    DETRAY_HOST_DEVICE
    constexpr bool is_on_portal() const { return cast_impl().on_portal(); }

    /// @returns @c true if the track has encountered a surface with material
    DETRAY_HOST_DEVICE
    constexpr bool encountered_sf_material() const {
        return is_on_surface() && (current().sf_desc.has_material());
    }

    /// @returns @c true if the navigation finished successfully
    DETRAY_HOST_DEVICE
    constexpr bool finished() const {
        return (status() == navigation::status::e_exit) && !is_alive();
    }

    /// @tparam surface_t the surface interface type that is required
    /// @returns current detector surface the navigator is on - const
    template <template <typename> class surface_t = tracking_surface>
    DETRAY_HOST_DEVICE constexpr auto current_surface() const {
        return surface_t{detector(), current().sf_desc};
    }

    /// @tparam volume_t the volume interface type that is required
    /// @returns current detector volume the navigation stream is in
    template <template <typename> class volume_t = tracking_volume>
    DETRAY_HOST_DEVICE constexpr auto current_volume() const {
        return volume_t{detector(), cast_impl().get_volume_idx()};
    }

    /// @tparam surface_t the surface interface type that is required
    /// @returns the next surface the navigator intends to reach
    template <template <typename> class surface_t = tracking_surface>
    DETRAY_HOST_DEVICE constexpr auto next_surface() const {
        return surface_t{detector(), target().sf_desc};
    }

    /// @tparam volume_t the volume interface type that is required
    /// @returns the next volume the navigator intends to reach
    template <template <typename> class volume_t = tracking_volume>
    DETRAY_HOST_DEVICE constexpr auto next_volume() const {
        return volume_t{detector(), target().volume_link};
    }

    /// @returns the navigation inspector - const
    DETRAY_HOST_DEVICE
    constexpr const auto& inspector() const {
        return cast_impl().get_inspector();
    }

    /// Navigation status that cannot be recovered from.
    DETRAY_HOST_DEVICE
    constexpr bool abort() {
        cast_impl().do_abort();
        // Call inspector. @Note This runs on the navigation state base type!
        cast_impl().get_inspector()(*this, "Aborted: ");

        return is_alive();
    }

    /// Navigation reaches final target or leaves detector world.
    DETRAY_HOST_DEVICE
    constexpr bool exit() {
        cast_impl().do_exit();
        // Call inspector. @Note This runs on the navigation state base type!
        cast_impl().get_inspector()(*this, "Exited: ");

        return is_alive();
    }
};

}  // namespace detray
