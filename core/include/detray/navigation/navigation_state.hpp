/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/detail/containers.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/barcode.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/utils/ranges.hpp"

namespace detray::navigation {

/// @enum NavigationDirection
/// The navigation direction is always with respect to a given track direction
enum class direction : std::int_least8_t { e_backward = -1, e_forward = 1 };

/// @enum Navigation status flags
enum class status : std::int_least8_t {
    e_abort = -3,          ///< error ocurred, navigation will be aborted
    e_exit = -2,           ///< navigation finished/reached the end of geometry
    e_unknown = -1,        ///< unknown state/not initialized
    e_towards_object = 0,  ///< move towards next geometry object
    e_on_object = 1,       ///< reached a geometry object that is not a portal
    e_on_portal = 2,       ///< reached portal (material) surface
};

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
/// The state is passed between navigation calls and is accessible to the
/// actors in the propagation, for which it defines the public interface
/// towards the navigation. The navigator class is responsible for updating the
/// elements in the navigation state with every navigation call, establishing
/// 'full trust' after changes to the track state reduced the trust level.
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
class state {

    protected:
    // Linear algebra types
    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using point3_t = dpoint3D<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    // Result of a geometry object evaluation
    using candidate_t = intersection_t;
    using candidate_cache_t = darray<candidate_t, k_cache_capacity>;
    using candidate_itr_t = typename candidate_cache_t::iterator;
    using candidate_const_itr_t = typename candidate_cache_t::const_iterator;
    using dist_t = std::int_least8_t;

    // Link type to the next detector volume
    using nav_link_type = typename detector_t::surface_type::navigation_link;

    public:
    using value_type = candidate_t;
    using detector_type = detector_t;

    // Vecmem memory management (in case an inspector needs dynamic memory)
    using view_type = detray::detail::get_view_t<inspector_t>;
    using const_view_type = detray::detail::get_view_t<const inspector_t>;

    /// Default constructor (needs a detector)
    state() = delete;

    /// Constructor using a given detector @param det
    DETRAY_HOST_DEVICE
    constexpr explicit state(const detector_type &det) : m_detector(&det) {}

    /// Construct from detector @param det and inspector view @param view
    template <concepts::device_view view_t>
    DETRAY_HOST_DEVICE constexpr state(const detector_type &det, view_t view)
        : m_detector(&det), m_inspector(view) {}

    /// @returns a reference to the detector
    DETRAY_HOST_DEVICE
    constexpr auto detector() const -> const detector_type & {
        return (*m_detector);
    }

    /// @returns the navigation heartbeat (is this stream still being updated?)
    DETRAY_HOST_DEVICE
    constexpr bool is_alive() const { return m_heartbeat; }

    /// @returns the navigation status (on object, towards object etc)
    DETRAY_HOST_DEVICE
    constexpr navigation::status status() { return m_status; }

    /// @returns the navigation trust level (how outdated is the state?)
    DETRAY_HOST_DEVICE
    constexpr navigation::trust_level trust_level() { return m_trust_level; }

    /// @returns currently cached candidates - const
    DETRAY_HOST_DEVICE
    constexpr auto candidates() const -> const candidate_cache_t & {
        return m_candidates;
    }

    /// @returns number of candidates the cache can hold - const
    DETRAY_HOST_DEVICE
    constexpr dindex n_candidates() const { return k_cache_capacity; }

    /// @returns current/previous object that was reached - const
    DETRAY_HOST_DEVICE
    constexpr auto current() const -> const candidate_t & {
        if constexpr (k_cache_capacity > 1) {
            assert(m_next > 0);
            return this->candidates()[static_cast<std::size_t>(m_next - 1)];
        } else {
            // Direct navigator
            return this->candidates()[static_cast<std::size_t>(m_next)];
        }
    }

    /// @returns next object that we want to reach (current target) - const
    DETRAY_HOST_DEVICE
    constexpr auto target() const -> const candidate_t & {
        return m_candidates[static_cast<std::size_t>(m_next)];
    }

    /// Scalar representation of the navigation state,
    /// @returns distance to next target
    DETRAY_HOST_DEVICE
    constexpr scalar_t operator()() const {
        return static_cast<scalar_t>(direction()) * target().path;
    }

    /// @returns current volume index - const
    DETRAY_HOST_DEVICE
    constexpr auto volume() const -> nav_link_type { return m_volume_index; }

    /// Set start/new volume from volume index @param v
    DETRAY_HOST_DEVICE
    constexpr void set_volume(dindex v) {
        assert(detail::is_invalid_value(static_cast<nav_link_type>(v)) ||
               v < detector().volumes().size());
        if (v != m_volume_index) {
            // Make sure the new volume is properly initialized
            set_no_trust();
        }
        m_volume_index = static_cast<nav_link_type>(v);
    }

    /// @returns barcode of the detector surface the navigator is on
    /// (invalid when not on surface) - const
    DETRAY_HOST_DEVICE
    constexpr auto barcode() const -> geometry::barcode {
        return current().sf_desc.barcode();
    }

    /// @returns current navigation status - const
    DETRAY_HOST_DEVICE
    constexpr auto status() const -> navigation::status { return m_status; }

    /// @returns current navigation direction - const
    DETRAY_HOST_DEVICE
    constexpr auto direction() const -> navigation::direction {
        return m_direction;
    }

    /// Set direction in which the navigator should search for candidates
    DETRAY_HOST_DEVICE
    constexpr void set_direction(const navigation::direction dir) {
        m_direction = dir;
    }

    /// @returns navigation trust level - const
    DETRAY_HOST_DEVICE
    constexpr auto trust_level() const -> navigation::trust_level {
        return m_trust_level;
    }

    /// Update navigation trust level to no trust
    DETRAY_HOST_DEVICE
    constexpr void set_no_trust() {
        m_trust_level = navigation::trust_level::e_no_trust;
    }

    /// Update navigation trust level to high trust
    DETRAY_HOST_DEVICE
    constexpr void set_high_trust() {
        m_trust_level = m_trust_level <= navigation::trust_level::e_high
                            ? m_trust_level
                            : navigation::trust_level::e_high;
    }

    /// Update navigation trust level to fair trust
    DETRAY_HOST_DEVICE
    constexpr void set_fair_trust() {
        m_trust_level = m_trust_level <= navigation::trust_level::e_fair
                            ? m_trust_level
                            : navigation::trust_level::e_fair;
    }

    /// @returns true if the track has reached a detector surface
    DETRAY_HOST_DEVICE
    constexpr bool is_on_surface() const {
        return (m_status >= navigation::status::e_on_object);
    }

    /// @returns true if the track has reached a sensitive surface
    DETRAY_HOST_DEVICE
    constexpr bool is_on_sensitive() const {
        return (m_status == navigation::status::e_on_object) &&
               current().sf_desc.is_sensitive();
    }

    /// @returns true if the track has reached a passive surface
    DETRAY_HOST_DEVICE
    constexpr bool is_on_passive() const {
        return (m_status == navigation::status::e_on_object) &&
               current().sf_desc.is_passive();
    }

    /// @returns true if the track has reached a portal surface
    DETRAY_HOST_DEVICE
    constexpr bool is_on_portal() const {
        return m_status == navigation::status::e_on_portal;
    }

    /// @returns true if the track has encountered a surface with material
    DETRAY_HOST_DEVICE
    constexpr bool encountered_sf_material() const {
        return is_on_surface() && (current().sf_desc.has_material());
    }

    /// @returns flag that indicates whether navigation was successful
    DETRAY_HOST_DEVICE
    constexpr bool is_complete() const {
        // Normal exit for this navigation?
        return (m_status == navigation::status::e_exit) && !m_heartbeat;
    }

    /// @param surface_t the surface interface type that is required
    /// @returns current detector surface the navigator is on - const
    template <template <typename> class surface_t = tracking_surface>
    DETRAY_HOST_DEVICE constexpr auto current_surface() const {
        assert(is_on_surface());
        return surface_t<detector_type>{*m_detector, current().sf_desc};
    }

    /// @param volume_t the volume interface type that is required
    /// @returns current detector volume of the navigation stream
    template <template <typename> class volume_t = tracking_volume>
    DETRAY_HOST_DEVICE constexpr auto current_volume() const {
        return volume_t<detector_type>{*m_detector, m_volume_index};
    }

    /// @param surface_t the surface interface type that is required
    /// @returns the next surface the navigator intends to reach
    template <template <typename> class surface_t = tracking_surface>
    DETRAY_HOST_DEVICE constexpr auto next_surface() const {
        return surface_t<detector_type>{*m_detector, target().sf_desc};
    }

    /// @param volume_t the volume interface type that is required
    /// @returns the next volume the navigator intends to reach
    template <template <typename> class volume_t = tracking_volume>
    DETRAY_HOST_DEVICE constexpr auto next_volume() const {
        return volume_t<detector_type>{*m_detector, target().sf_desc.volume()};
    }

    /// @returns the navigation inspector - const
    DETRAY_HOST_DEVICE
    constexpr const auto &inspector() const { return m_inspector; }

    /// @returns the navigation inspector
    DETRAY_HOST_DEVICE
    constexpr auto &inspector() { return m_inspector; }

    protected:
    /// Set the heartbeat to @param b
    DETRAY_HOST_DEVICE
    constexpr void heartbeat(bool b) { m_heartbeat = b; }

    /// Set the status to @param s
    DETRAY_HOST_DEVICE
    constexpr void status(navigation::status s) { m_status = s; }

    /// Reset the trustlevel to @param t
    DETRAY_HOST_DEVICE
    constexpr void trust_level(navigation::trust_level t) { m_trust_level = t; }

    /// @returns the index to the target surface
    DETRAY_HOST_DEVICE
    constexpr dist_t next_index() const { return m_next; }

    /// Set the next surface that we want to reach (update target)
    DETRAY_HOST_DEVICE
    constexpr void next_index(dindex pos) {
        m_next = static_cast<dist_t>(pos);
        assert(m_next < static_cast<dist_t>(k_cache_capacity));
    }

    /// Set the next surface that we want to reach (update target)
    DETRAY_HOST_DEVICE
    constexpr auto advance() -> void {
        ++this->m_next;
        assert(this->m_next < static_cast<dist_t>(k_cache_capacity) + 1);
    }

    /// @returns currently cached candidates
    DETRAY_HOST_DEVICE
    constexpr auto candidates() -> candidate_cache_t & { return m_candidates; }

    /// @returns current/previous object that was reached
    DETRAY_HOST_DEVICE
    constexpr auto current() -> candidate_t & {
        if constexpr (k_cache_capacity > 1) {
            assert(m_next > 0);
            return this->candidates()[static_cast<std::size_t>(m_next - 1)];
        } else {
            // Direct navigator
            return this->candidates()[static_cast<std::size_t>(m_next)];
        }
    }

    /// @returns next object that we want to reach (current target)
    DETRAY_HOST_DEVICE
    constexpr auto target() -> candidate_t & {
        assert(static_cast<std::size_t>(m_next) < m_candidates.size());
        return m_candidates[static_cast<std::size_t>(m_next)];
    }

    /// @returns true if a candidate lies on a surface - const
    DETRAY_HOST_DEVICE constexpr bool has_reached_surface(
        const candidate_t &candidate, const navigation::config &cfg) const {
        return (math::fabs(candidate.path) < cfg.path_tolerance);
    }

    /// Clear the state
    DETRAY_HOST_DEVICE
    constexpr void clear() {
        // Mark all data in the cache as unreachable
        for (std::size_t i = 0u; i < k_cache_capacity; ++i) {
            m_candidates[i].path = std::numeric_limits<scalar_t>::max();
        }
        m_next = 0;
    }

    /// Navigation state that cannot be recovered from. Leave the state
    /// data for inspection.
    DETRAY_HOST_DEVICE
    constexpr void abort() {
        m_status = navigation::status::e_abort;
        m_heartbeat = false;
        // Don't do anything if aborted
        m_trust_level = navigation::trust_level::e_full;
    }

    /// Navigation reaches final target or leaves detector world. Stop
    /// navigation.
    DETRAY_HOST_DEVICE
    constexpr void exit() {
        m_status = navigation::status::e_exit;
        m_heartbeat = false;
        m_trust_level = navigation::trust_level::e_full;
    }

    private:
    /// Our cache of candidates (intersections with any kind of surface)
    candidate_cache_t m_candidates;

    /// Detector pointer
    const detector_type *m_detector{nullptr};

    /// Index of current navigation volume (in the detector volume container)
    nav_link_type m_volume_index{0u};

    /// The next best candidate (target) in the candidate cache
    dist_t m_next{0};

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
