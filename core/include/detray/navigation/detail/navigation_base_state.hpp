/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/navigation/navigation_state_interface.hpp"

namespace detray::navigation {

/// @brief A void inpector that does nothing.
///
/// Inspectors can be plugged in to understand the current navigation state.
struct void_inspector {

    struct void_view : public detray::detail::dbase_view {};

    using view_type = void_view;
    using const_view_type = const void_view;

    constexpr void_inspector() = default;

    DETRAY_HOST_DEVICE
    constexpr explicit void_inspector(
        const void_view & /*ignored*/) { /*Do nothing*/
    }

    template <typename state_t>
    DETRAY_HOST_DEVICE constexpr void operator()(
        const state_t & /*ignored*/, const char * /*ignored*/) const {
        /*Do nothing*/
    }
};

/// @brief A navigation state object used to cache the information of the
/// current navigation stream.
///
/// This is the base class for all concrete navigation streams that are defined
/// in the navigator implementations. It inherits its public interface from the
/// @c navigation_state_interface class.
/// The state is passed between navigation calls and is updated by its
/// navigator, establishing 'full trust' after changes to the track state
/// may have reduced the trust level.
///
/// @tparam detector_t the type of the detector that is being navigated
/// @tparam k_cache_capacity number of object candidates the state can hold
/// @tparam inspector_t is a validation inspector that can record information
///         about the navigation state at different points of the nav. flow.
/// @tparam intersection_t result of an intersection operation
template <typename detector_t, std::size_t k_cache_capacity,
          typename inspector_t = navigation::void_inspector,
          typename intersection_t =
              intersection2D<typename detector_t::surface_type,
                             typename detector_t::algebra_type, false>>
class base_state
    : public navigation_state_interface<base_state<
          detector_t, k_cache_capacity, inspector_t, intersection_t>> {

    friend class navigation_state_interface<base_state>;

    using public_interface_t = navigation_state_interface<base_state>;

    protected:
    // Linear algebra types
    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;

    // Result of a geometry object evaluation
    using candidate_t = intersection_t;
    using candidate_cache_t = darray<candidate_t, k_cache_capacity>;
    using candidate_itr_t = typename candidate_cache_t::iterator;
    using candidate_const_itr_t = typename candidate_cache_t::const_iterator;
    using dist_t = std::int_least8_t;

    public:
    // Link type to the next detector volume
    using nav_link_type = typename detector_t::surface_type::navigation_link;

    using value_type = candidate_t;
    using detector_type = detector_t;
    using inspector_type = inspector_t;

    // Vecmem memory management (in case an inspector needs dynamic memory)
    using view_type = detray::detail::get_view_t<inspector_t>;
    using const_view_type = detray::detail::get_view_t<const inspector_t>;

    /// No default constructor (needs a detector)
    base_state() = delete;

    /// Constructor using a given detector @param det
    DETRAY_HOST_DEVICE
    constexpr explicit base_state(const detector_type &det)
        : m_detector(&det) {}

    /// Construct from detector @param det and inspector view @param view
    template <concepts::device_view view_t>
    DETRAY_HOST_DEVICE constexpr base_state(const detector_type &det,
                                            view_t view)
        : m_detector(&det), m_inspector(view) {}

    /// Make the constant range access from the public interface accessible to
    /// overload resolution
    /// @{
    using public_interface_t::begin;
    using public_interface_t::end;
    /// @}

    protected:
    /// @returns a reference to the detector that is being navigated - const
    DETRAY_HOST_DEVICE
    constexpr auto get_detector() const -> const detector_type & {
        return *m_detector;
    }

    /// @returns the navigation heartbeat (is this stream still being updated?)
    DETRAY_HOST_DEVICE
    constexpr bool get_heartbeat() const { return m_heartbeat; }

    /// Sets the navigation heartbeat to @param h
    DETRAY_HOST_DEVICE
    constexpr void set_heartbeat(bool h) { m_heartbeat = h; }

    /// @returns the navigation status (on object, towards object etc)
    DETRAY_HOST_DEVICE
    constexpr auto get_status() const -> navigation::status { return m_status; }

    /// Sets the navigation status to @param s
    DETRAY_HOST_DEVICE
    constexpr void set_status(navigation::status s) { m_status = s; }

    /// @returns the navigation trust level (how outdated is the state?)
    DETRAY_HOST_DEVICE
    constexpr auto get_trust_level() const -> navigation::trust_level {
        return m_trust_level;
    }

    /// Sets the navigation trust level to @param level
    DETRAY_HOST_DEVICE
    constexpr void set_trust_level(navigation::trust_level level) {
        m_trust_level = level;
    }

    /// @returns current navigation direction - const
    DETRAY_HOST_DEVICE
    constexpr auto get_direction() const -> navigation::direction {
        return m_direction;
    }

    /// Change direction in which the navigator should search for candidates
    DETRAY_HOST_DEVICE
    constexpr void set_direction(const navigation::direction dir) {
        m_direction = dir;
    }

    /// @returns start position of valid candidate range - const
    DETRAY_HOST_DEVICE
    constexpr auto get_begin() const -> candidate_const_itr_t {
        candidate_const_itr_t itr = m_candidates.begin();

        if constexpr (k_cache_capacity > 1) {
            detray::ranges::advance(
                itr, (on_surface() && (m_next >= 1)) ? m_next - 1 : m_next);
        }

        return itr;
    }

    /// @returns start position of valid candidate range - non-const
    DETRAY_HOST_DEVICE
    constexpr auto get_begin() -> candidate_itr_t {
        candidate_itr_t itr = m_candidates.begin();

        if constexpr (k_cache_capacity > 1) {
            detray::ranges::advance(
                itr, (on_surface() && (m_next >= 1)) ? m_next - 1 : m_next);
        }

        return itr;
    }

    /// @returns sentinel of the valid candidate range - const
    DETRAY_HOST_DEVICE
    constexpr auto get_end() const -> candidate_const_itr_t {
        candidate_const_itr_t itr = m_candidates.begin();
        detray::ranges::advance(itr, m_last + 1);
        return itr;
    }

    /// @returns sentinel of the valid candidate range - non-const
    DETRAY_HOST_DEVICE
    constexpr auto get_end() -> candidate_itr_t {
        candidate_itr_t itr = m_candidates.begin();
        detray::ranges::advance(itr, m_last + 1);
        return itr;
    }

    /// Non-const range-access (only available from within the derived states)
    /// @{
    DETRAY_HOST_DEVICE
    constexpr auto begin() { return get_begin(); }
    DETRAY_HOST_DEVICE
    constexpr auto end() { return get_end(); }
    /// @}

    /// @returns current navigation direction - const
    DETRAY_HOST_DEVICE
    constexpr dist_t next_index() const { return m_next; }

    /// Set the next surface that we want to reach (update target)
    DETRAY_HOST_DEVICE
    constexpr void set_next_index(dist_t n) {
        assert(n < static_cast<dist_t>(k_cache_capacity));
        m_next = n;
    }

    /// @returns current navigation direction - const
    DETRAY_HOST_DEVICE
    constexpr dist_t last_index() const { return m_last; }

    /// Updates the position of the last valid candidate
    DETRAY_HOST_DEVICE
    constexpr void set_last_index(dist_t l) {
        assert(m_next <= l + 1);
        assert(l < static_cast<dist_t>(k_cache_capacity));
        m_last = l;
    }

    /// @returns currently cached candidates - const
    DETRAY_HOST_DEVICE
    constexpr const candidate_cache_t &get_candidates() const {
        return m_candidates;
    }

    /// @returns currently cached candidates - const
    DETRAY_HOST_DEVICE
    constexpr candidate_cache_t &get_candidates() { return m_candidates; }

    /// @returns number of candidates the cache can hold - const
    DETRAY_HOST_DEVICE
    consteval dindex get_capacity() const { return k_cache_capacity; }

    /// @returns number of candidates the cache currently holds - const
    DETRAY_HOST_DEVICE
    dindex get_n_candidates() const {
        assert(m_last - m_next + 1 >= 0);
        return static_cast<dindex>(m_last - m_next + 1);
    }

    /// @returns current volume index - const
    DETRAY_HOST_DEVICE
    constexpr nav_link_type get_volume_idx() const { return m_volume_index; }

    /// Set start/new volume from volume index @param v
    DETRAY_HOST_DEVICE
    constexpr void set_volume_idx(dindex v) {
        assert(detail::is_invalid_value(static_cast<nav_link_type>(v)) ||
               v < m_detector->volumes().size());

        const auto v_link{static_cast<nav_link_type>(v)};
        if (v_link != m_volume_index) {
            // Make sure the new volume is properly initialized
            m_trust_level = navigation::trust_level::e_no_trust;
        }
        m_volume_index = v_link;
    }

    /// @returns current candidate that was reached - const
    DETRAY_HOST_DEVICE
    constexpr auto get_current() const -> const candidate_t & {
        if constexpr (k_cache_capacity > 1) {
            assert(m_next > 0);
            return m_candidates[static_cast<std::size_t>(m_next - 1)];
        } else {
            return m_candidates[0];
        }
    }

    /// @returns next candidate that we want to reach - const
    DETRAY_HOST_DEVICE
    constexpr auto get_target() const -> const candidate_t & {
        return m_candidates[static_cast<std::size_t>(m_next)];
    }

    /// @returns next candidate that we want to reach - non- const
    DETRAY_HOST_DEVICE
    constexpr auto get_target() -> candidate_t & {
        return m_candidates[static_cast<std::size_t>(m_next)];
    }

    /// @returns @c true if the track has reached a detector surface
    DETRAY_HOST_DEVICE
    constexpr bool on_surface() const {
        return (m_status >= navigation::status::e_on_object);
    }

    /// @returns @c true if the track has reached a portal surface
    DETRAY_HOST_DEVICE
    constexpr bool on_portal() const {
        return m_status == navigation::status::e_on_portal;
    }

    /// @returns the navigation inspector - const
    DETRAY_HOST_DEVICE
    constexpr const auto &get_inspector() const { return m_inspector; }

    /// @returns the navigation inspector - non-const
    DETRAY_HOST_DEVICE
    constexpr auto &get_inspector() { return m_inspector; }

    /// Clear the state
    DETRAY_HOST_DEVICE
    constexpr void clear_cache() {
        // Mark all data in the cache as unreachable
        for (std::size_t i = 0u; i < k_cache_capacity; ++i) {
            m_candidates[i].path = std::numeric_limits<scalar_t>::max();
        }
        m_next = 0;
        m_last = -1;
    }

    /// Navigation state that cannot be recovered from. Leave the state
    /// data for inspection.
    DETRAY_HOST_DEVICE
    constexpr void do_abort() {
        m_status = navigation::status::e_abort;
        m_heartbeat = false;
        // Don't do anything if aborted
        m_trust_level = navigation::trust_level::e_full;

        clear_cache();
    }

    /// Navigation reaches final target or leaves detector world. Stop
    /// navigation.
    DETRAY_HOST_DEVICE
    constexpr void do_exit() {
        m_status = navigation::status::e_exit;
        m_heartbeat = false;
        m_trust_level = navigation::trust_level::e_full;
    }

    /// Our cache of candidates (intersections with any kind of surface)
    candidate_cache_t m_candidates;

    /// Detector pointer
    const detector_type *m_detector{nullptr};

    /// Index of current navigation volume (in the detector volume container)
    nav_link_type m_volume_index{0u};

    /// The next best candidate (target) in the candidate cache
    dist_t m_next{0};

    /// The last reachable candidate: m_last < k_cache_capacity
    /// Can never be advanced beyond the last element
    dist_t m_last{-1};

    /// The navigation status
    navigation::status m_status{navigation::status::e_unknown};

    /// The navigation trust level determines how this states cache is to
    /// be updated in the current navigation call
    navigation::trust_level m_trust_level{navigation::trust_level::e_no_trust};

    /// The navigation direction
    navigation::direction m_direction{navigation::direction::e_forward};

    /// Heartbeat of this navigation flow signals navigation is alive
    bool m_heartbeat{false};

    /// The inspector type of this navigation engine
    [[no_unique_address]] inspector_t m_inspector;
};

}  // namespace detray::navigation
