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

namespace navigation {

static constexpr std::size_t default_sequence_size{100u};

}

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

    class state : public detray::ranges::view_interface<state> {

        using candidate_t = intersection_type;

        public:
        using sequence_t =
            vecmem::device_vector<detray::geometry::barcode::value_t>;
        using detector_type = direct_navigator::detector_type;
        using nav_link_type =
            typename detector_type::surface_type::navigation_link;
        using view_type = bool;

        state() = delete;

        template <typename vector_t>
        DETRAY_HOST_DEVICE explicit state(const detector_t &det,
                                          const vector_t &sequence)
            : m_detector(&det), m_sequence(sequence) {

            // TODO: Copy the sequence
            m_it = m_sequence.begin();
            m_it_rev = m_sequence.rbegin();
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

        /// @returns next object that we want to reach (current target) - const
        DETRAY_HOST_DEVICE
        inline auto target() -> candidate_t & { return m_candidate; }

        /// @returns next object that we want to reach (current target) - const
        DETRAY_HOST_DEVICE
        void update() {
            m_candidate.sf_desc.set_barcode(get_target_barcode());
            set_volume(m_candidate.volume_link);
            if (!is_init()) {
                m_candidate_prev.sf_desc.set_barcode(get_current_barcode());
            }
        }

        /// @returns current detector surface the navigator is on
        /// (cannot be used when not on surface) - const
        DETRAY_HOST_DEVICE
        inline auto get_surface() const {
            assert(is_on_surface());
            return tracking_surface<detector_type>{*m_detector,
                                                   current().sf_desc};
        }
        /*
        DETRAY_HOST_DEVICE
        inline auto target() const -> const candidate_t & {
            return m_candidate;
        }
        */

        /// @returns current navigation status - const
        DETRAY_HOST_DEVICE
        inline auto status() const -> navigation::status { return m_status; }

        DETRAY_HOST_DEVICE
        inline bool is_init() {
            if (m_direction == navigation::direction::e_forward) {
                return m_it == m_sequence.begin();
            } else {
                return m_it_rev == m_sequence.rbegin();
            }
        }

        /// @return the reference of track state pointed by the iterator
        DETRAY_HOST_DEVICE
        const detray::geometry::barcode get_target_barcode() const {
            if (m_direction == navigation::direction::e_forward) {
                return detray::geometry::barcode{*m_it};
            } else {
                return detray::geometry::barcode{*m_it_rev};
            }
        }

        /// @return the reference of track state pointed by the iterator
        DETRAY_HOST_DEVICE
        const detray::geometry::barcode get_current_barcode() const {
            if (m_direction == navigation::direction::e_forward) {
                return detray::geometry::barcode{*(m_it - 1)};
            } else {
                return detray::geometry::barcode{*(m_it_rev - 1)};
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
        bool is_complete() {
            if ((m_direction == navigation::direction::e_forward) &&
                m_it == m_sequence.end()) {
                return true;
            } else if ((m_direction != navigation::direction::e_forward) &&
                       m_it_rev == m_sequence.rend()) {
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
            return (is_on_surface()) && (target().sf_desc.material().id() !=
                                         detector_t::materials::id::e_none);
        }

        /// Helper method to check the track has reached a sensitive surface
        DETRAY_HOST_DEVICE
        inline auto is_on_sensitive() const -> bool {
            return (m_status == navigation::status::e_on_module) &&
                   (barcode().id() == surface_id::e_sensitive);
        }

        DETRAY_HOST_DEVICE
        inline auto barcode() const -> geometry::barcode {
            return m_candidate.sf_desc.barcode();
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

        /// @returns current detector volume of the navigation stream
        DETRAY_HOST_DEVICE
        inline auto get_volume() const {
            return tracking_volume<detector_type>{*m_detector, m_volume_index};
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
        /// @TODO: Make it array with multiple candidates, and resort the target
        /// elements based on the path length
        candidate_t m_candidate;
        candidate_t m_candidate_prev;

        /// Detector pointer
        const detector_type *m_detector{nullptr};

        /// Index in the detector volume container of current navigation volume
        nav_link_type m_volume_index{0u};

        /// Target surfaces
        sequence_t m_sequence;

        // iterator for forward direction
        typename sequence_t::iterator m_it;

        // iterator for backward direction
        typename sequence_t::reverse_iterator m_it_rev;

        /// The navigation direction
        navigation::direction m_direction{navigation::direction::e_forward};

        /// The navigation status
        navigation::status m_status{navigation::status::e_unknown};

        /// Heartbeat of this navigation flow signals navigation is alive
        bool m_heartbeat{false};
    };

    template <typename track_t>
    DETRAY_HOST_DEVICE inline void init(const track_t &track, state &navigation,
                                        const navigation::config &cfg,
                                        const context_type &ctx) const {
        navigation.m_heartbeat = true;

        // Set the geometry barcode for the candidate
        navigation.update();

        return;
    }

    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool update(const track_t &track,
                                          state &navigation,
                                          const navigation::config &cfg,
                                          const context_type &ctx = {}) const {

        if (navigation.get_target_barcode().is_invalid()) {
            return false;
        }

        if (navigation.is_on_surface(navigation.target(), cfg)) {
            navigation.next();
            navigation.update();
            navigation.m_status = (navigation.target().sf_desc.is_portal())
                                      ? navigation::status::e_on_portal
                                      : navigation::status::e_on_module;
        } else {
            // Otherwise the track is moving towards a surface
            navigation.m_status = navigation::status::e_towards_object;
        }

        const auto &det = navigation.detector();
        const auto sf = tracking_surface{det, navigation.target().sf_desc};

        // Check whether this candidate is reachable by the track
        return sf.template visit_mask<intersection_update<ray_intersector>>(
            detail::ray<algebra_type>(
                track.pos(),
                static_cast<scalar_type>(navigation.direction()) * track.dir()),
            navigation.target(), det.transform_store(), ctx,
            sf.is_portal() ? darray<scalar_type, 2>{0.f, 0.f}
                           : darray<scalar_type, 2>{cfg.min_mask_tolerance,
                                                    cfg.max_mask_tolerance},
            static_cast<scalar_type>(cfg.mask_tolerance_scalor),
            static_cast<scalar_type>(cfg.overstep_tolerance));
    }
};

}  // namespace detray