/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/algorithms.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/barcode.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/intersection_kernel.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/ranges.hpp"

namespace detray {

template <typename detector_t>
class direct_navigator {

    public:
    using detector_type = detector_t;
    using context_type = detector_type::geometry_context;
    using algebra_type = typename detector_type::algebra_type;
    using scalar_type = dscalar<algebra_type>;
    using intersection_type =
        intersection2D<typename detector_t::surface_type,
                       typename detector_t::algebra_type, false>;

    class state {

        friend struct intersection_update<ray_intersector>;

        using candidate_t = intersection_type;

        public:
        using sequence_t = vecmem::device_vector<detray::geometry::barcode>;
        using detector_type = direct_navigator::detector_type;
        using nav_link_type =
            typename detector_type::surface_type::navigation_link;
        using view_type = dvector_view<detray::geometry::barcode>;

        state() = delete;

        DETRAY_HOST_DEVICE explicit state(const detector_t &det,
                                          const view_type &sequence)
            : m_detector(&det), m_sequence(sequence) {

            m_it = m_sequence.cbegin();
            m_it_rev = m_sequence.crbegin();
        }

        /// Scalar representation of the navigation state,
        /// @returns distance to next
        DETRAY_HOST_DEVICE
        scalar_type operator()() const {
            return static_cast<scalar_type>(direction()) * target().path;
        }

        /// @returns a pointer of detector
        DETRAY_HOST_DEVICE
        const detector_type &detector() const { return (*m_detector); }

        /// @returns the navigation heartbeat
        DETRAY_HOST_DEVICE
        bool is_alive() const { return m_heartbeat; }

        /// @returns current/previous object that was reached
        DETRAY_HOST_DEVICE
        inline auto current() const -> const candidate_t & {
            return m_candidate_prev;
        }

        /// @returns next object that we want to reach (current target) - const
        DETRAY_HOST_DEVICE
        inline auto target() const -> const candidate_t & {
            return m_candidate;
        }

        DETRAY_HOST_DEVICE
        inline auto target() -> candidate_t & { return m_candidate; }

        DETRAY_HOST_DEVICE
        void update_candidate(bool update_candidate_prev = true) {

            if (update_candidate_prev) {
                m_candidate_prev = m_candidate;
            }

            if (!is_complete()) {
                m_candidate.sf_desc = m_detector->surface(get_target_barcode());
                m_candidate.volume_link =
                    tracking_surface{*m_detector, m_candidate.sf_desc}
                        .volume_link();
                m_candidate.path = std::numeric_limits<scalar_type>::max();
                set_volume(m_candidate.volume_link);
            }
        }

        /// @returns current detector surface the navigator is on
        /// (cannot be used when not on surface) - const
        DETRAY_HOST_DEVICE
        inline auto get_surface() const -> tracking_surface<detector_type> {
            assert(is_on_surface());
            return tracking_surface<detector_type>{*m_detector,
                                                   current().sf_desc};
        }

        /// @returns current navigation status - const
        DETRAY_HOST_DEVICE
        inline auto status() const -> navigation::status { return m_status; }

        DETRAY_HOST_DEVICE
        inline bool is_init() const {
            if (m_direction == navigation::direction::e_forward) {
                return m_it == m_sequence.cbegin();
            } else {
                return m_it_rev == m_sequence.crbegin();
            }
        }

        DETRAY_HOST_DEVICE
        detray::geometry::barcode get_target_barcode() const {
            if (m_direction == navigation::direction::e_forward) {
                return *m_it;
            } else {
                return *m_it_rev;
            }
        }

        DETRAY_HOST_DEVICE
        detray::geometry::barcode get_current_barcode() const {
            if (m_direction == navigation::direction::e_forward) {
                return *(m_it - 1);
            } else {
                return *(m_it_rev - 1);
            }
        }

        /// Advance the iterator
        DETRAY_HOST_DEVICE
        void next() {
            if (m_direction == navigation::direction::e_forward) {
                m_it++;
            } else {
                m_it_rev++;
            }
        }

        /// @return true if the iterator reaches the end of vector
        DETRAY_HOST_DEVICE
        bool is_complete() const {
            if ((m_direction == navigation::direction::e_forward) &&
                m_it == m_sequence.cend()) {
                return true;
            } else if ((m_direction == navigation::direction::e_backward) &&
                       m_it_rev == m_sequence.crend()) {
                return true;
            }
            return false;
        }

        /// @returns current navigation direction - const
        DETRAY_HOST_DEVICE
        inline auto direction() const -> navigation::direction {
            return m_direction;
        }

        /// Helper method to check the track has reached a module surface
        DETRAY_HOST_DEVICE
        inline auto is_on_surface() const -> bool {
            return (m_status == navigation::status::e_on_module ||
                    m_status == navigation::status::e_on_portal);
        }

        /// Helper method to check if a candidate lies on a surface - const
        DETRAY_HOST_DEVICE inline auto is_on_surface(
            const intersection_type &candidate,
            const navigation::config &cfg) const -> bool {
            return (math::fabs(candidate.path) < cfg.path_tolerance);
        }

        /// Helper method to check the track has encountered material
        DETRAY_HOST_DEVICE
        inline auto encountered_sf_material() const -> bool {
            return (is_on_surface()) && (current().sf_desc.material().id() !=
                                         detector_t::materials::id::e_none);
        }

        /// Helper method to check the track has reached a sensitive surface
        DETRAY_HOST_DEVICE
        inline auto is_on_sensitive() const -> bool {
            return (m_status == navigation::status::e_on_module) &&
                   (get_current_barcode().id() == surface_id::e_sensitive);
        }

        DETRAY_HOST_DEVICE
        inline auto barcode() const -> geometry::barcode {
            return m_candidate_prev.sf_desc.barcode();
        }

        /// @returns current volume (index) - const
        DETRAY_HOST_DEVICE
        inline auto volume() const -> nav_link_type { return m_volume_index; }

        /// Set start/new volume
        DETRAY_HOST_DEVICE
        inline void set_volume(dindex v) {
            assert(detail::is_invalid_value(static_cast<nav_link_type>(v)) ||
                   v < detector().volumes().size());
            m_volume_index = static_cast<nav_link_type>(v);
        }

        DETRAY_HOST_DEVICE
        inline auto abort(const char * = nullptr) -> bool {
            m_status = navigation::status::e_abort;
            m_heartbeat = false;
            return m_heartbeat;
        }

        template <typename debug_msg_generator_t>
        DETRAY_HOST_DEVICE inline auto abort(const debug_msg_generator_t &)
            -> bool {
            return abort();
        }

        DETRAY_HOST_DEVICE
        inline auto exit() -> bool {
            m_status = navigation::status::e_on_target;
            m_heartbeat = false;
            return m_heartbeat;
        }

        DETRAY_HOST_DEVICE
        inline auto pause() const -> bool { return false; }

        /// @returns current detector volume of the navigation stream
        DETRAY_HOST_DEVICE
        inline auto get_volume() const {
            return tracking_volume<detector_type>{*m_detector, m_volume_index};
        }

        /// Set direction
        DETRAY_HOST_DEVICE
        inline void set_direction(const navigation::direction dir) {
            m_direction = dir;
        }

        DETRAY_HOST_DEVICE
        inline void set_no_trust() { return; }

        DETRAY_HOST_DEVICE
        inline void set_full_trust() { return; }

        DETRAY_HOST_DEVICE
        inline void set_high_trust() { return; }

        DETRAY_HOST_DEVICE
        inline void set_fair_trust() { return; }

        /// Intersection candidate
        candidate_t m_candidate;
        candidate_t m_candidate_prev;

        /// Detector pointer
        const detector_type *m_detector{nullptr};

        /// Index in the detector volume container of current navigation volume
        nav_link_type m_volume_index{0u};

        /// Target surfaces
        sequence_t m_sequence;

        // iterator for forward direction
        typename sequence_t::const_iterator m_it;

        // iterator for backward direction
        typename sequence_t::const_reverse_iterator m_it_rev;

        /// The navigation direction
        navigation::direction m_direction{navigation::direction::e_forward};

        /// The navigation status
        navigation::status m_status{navigation::status::e_unknown};

        /// Step size when the valid intersection is not found for the target
        scalar_type safe_step_size = 10.f * unit<scalar_type>::mm;

        /// Heartbeat of this navigation flow signals navigation is alive
        bool m_heartbeat{false};
    };

    template <typename track_t>
    DETRAY_HOST_DEVICE inline void init(const track_t &track, state &navigation,
                                        const navigation::config &cfg,
                                        const context_type &ctx) const {
        // Do not resurrect a failed/finished navigation state
        assert(navigation.status() > navigation::status::e_on_target);

        if (navigation.is_complete()) {
            navigation.m_heartbeat = false;
            return;
        }

        navigation.m_heartbeat = true;
        navigation.update_candidate(!navigation.is_init());
        update(track, navigation, cfg, ctx);
    }

    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool update(
        const track_t &track, state &navigation, const navigation::config &cfg,
        const context_type &ctx = {},
        const bool is_before_actor_run = true) const {

        if (navigation.is_complete()) {
            navigation.m_heartbeat = false;
            return true;
        }

        assert(!navigation.get_target_barcode().is_invalid());
        update_intersection(track, navigation, cfg, ctx);

        if (is_before_actor_run) {
            if (navigation.is_on_surface(navigation.target(), cfg)) {
                navigation.m_status = (navigation.target().sf_desc.is_portal())
                                          ? navigation::status::e_on_portal
                                          : navigation::status::e_on_module;
                navigation.next();
                navigation.update_candidate(true);
                assert(navigation.is_on_surface(navigation.current(), cfg));

                if (!navigation.is_complete()) {
                    update_intersection(track, navigation, cfg, ctx);
                }

                // Return true to reset the step size
                return true;
            }
        }

        // Otherwise the track is moving towards a surface
        navigation.m_status = navigation::status::e_towards_object;

        // Return false to scale the step size with RK4
        return false;
    }

    template <typename track_t>
    DETRAY_HOST_DEVICE inline void update_intersection(
        const track_t &track, state &navigation, const navigation::config &cfg,
        const context_type &ctx = {}) const {

        const auto &det = navigation.detector();
        const auto sf = tracking_surface{det, navigation.target().sf_desc};

        const bool res =
            sf.template visit_mask<intersection_update<ray_intersector>>(
                detail::ray<algebra_type>(
                    track.pos(),
                    static_cast<scalar_type>(navigation.direction()) *
                        track.dir()),
                navigation.target(), det.transform_store(), ctx,
                darray<scalar_type, 2>{cfg.min_mask_tolerance,
                                       cfg.max_mask_tolerance},
                static_cast<scalar_type>(cfg.mask_tolerance_scalor),
                static_cast<scalar_type>(cfg.overstep_tolerance));

        // If an intersection is not found, proceed the track with safe step
        // size
        if (!res) {
            const auto path = navigation.target().path;
            navigation.update_candidate(false);
            navigation.target().path =
                math::copysign(navigation.safe_step_size, path);
        }
    }
};

}  // namespace detray
