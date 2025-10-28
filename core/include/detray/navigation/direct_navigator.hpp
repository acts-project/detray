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
#include "detray/definitions/navigation.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/barcode.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/intersection_kernel.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/navigation/navigation_state.hpp"
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
    using intersection_type = intersection2D<typename detector_t::surface_type,
                                             typename detector_t::algebra_type,
                                             !intersection::contains_pos>;

    class state : public navigation::base_state<state, detector_type, 2u,
                                                navigation::void_inspector,
                                                intersection_type> {
        friend class direct_navigator;
        friend struct intersection_update<ray_intersector>;

        using base_type = navigation::base_state<state, detector_type, 2u,
                                                 navigation::void_inspector,
                                                 intersection_type>;

        using candidate_t = intersection_type;
        using nav_link_t =
            typename detector_type::surface_type::navigation_link;

        public:
        using value_type = candidate_t;
        using sequence_type = vecmem::device_vector<detray::geometry::barcode>;

        using view_type = dvector_view<detray::geometry::barcode>;
        using const_view_type = dvector_view<const detray::geometry::barcode>;

        state() = delete;

        DETRAY_HOST_DEVICE constexpr state(const detector_t &det,
                                           const view_type &sequence)
            : base_type(det), m_sequence(sequence) {

            m_it = m_sequence.cbegin();
            m_it_rev = m_sequence.crbegin();

            // The next candidate is always stored in the second cache entry
            this->next_index(1u);
            this->last_index(1u);
        }

        /// @returns the direct navigator always has only one condidate
        DETRAY_HOST_DEVICE
        constexpr auto n_candidates() const -> dindex { return 1u; }

        DETRAY_HOST_DEVICE
        constexpr void update_candidate(bool update_candidate_prev = true) {
            if (update_candidate_prev) {
                this->candidates()[0] = this->candidates()[1];
            }

            if (!no_next_external()) {
                this->candidates()[1].sf_desc =
                    this->detector().surface(get_target_barcode());
                this->candidates()[1].volume_link =
                    detail::invalid_value<nav_link_t>();
                this->candidates()[1].set_path(
                    std::numeric_limits<scalar_type>::max());
                this->set_volume(this->candidates()[1].sf_desc.volume());
            }
        }

        DETRAY_HOST_DEVICE
        constexpr bool is_init() const {
            return is_forward() ? m_it == m_sequence.cbegin()
                                : m_it_rev == m_sequence.crbegin();
        }

        DETRAY_HOST_DEVICE
        constexpr detray::geometry::barcode get_target_barcode() const {
            return is_forward() ? *m_it : *m_it_rev;
        }

        DETRAY_HOST_DEVICE
        constexpr detray::geometry::barcode get_current_barcode() const {
            return is_forward() ? *(m_it - 1) : *(m_it_rev - 1);
        }

        /// Advance the iterator
        DETRAY_HOST_DEVICE
        constexpr void next_external() {
            if (is_forward()) {
                m_it++;
            } else {
                m_it_rev++;
            }
        }

        /// @return true if the iterator reaches the end of vector
        DETRAY_HOST_DEVICE
        constexpr bool no_next_external() const {
            return (is_forward() && m_it == m_sequence.cend()) ||
                   (!is_forward() && m_it_rev == m_sequence.crend());
        }

        DETRAY_HOST_DEVICE
        constexpr scalar_type safe_step_size() const {
            assert(m_safe_step_size > 0.f);
            return m_safe_step_size;
        }

        private:
        /// @return 'true' if the navigation direction is 'forward'
        DETRAY_HOST_DEVICE
        constexpr bool is_forward() const {
            return this->direction() == navigation::direction::e_forward;
        }

        /// Target surfaces
        sequence_type m_sequence;

        // iterator for forward direction
        typename sequence_type::const_iterator m_it;

        // iterator for backward direction
        typename sequence_type::const_reverse_iterator m_it_rev;

        /// Step size when the valid intersection is not found for the target
        scalar_type m_safe_step_size = 10.f * unit<scalar_type>::mm;
    };

    template <typename track_t>
    DETRAY_HOST_DEVICE inline void init(const track_t &track, state &navigation,
                                        const navigation::config &cfg,
                                        const context_type &ctx) const {
        DETRAY_VERBOSE_HOST_DEVICE("Called 'init()'");

        // Do not resurrect a failed/finished navigation state
        assert(navigation.is_alive());
        assert(!track.is_invalid());

        if (navigation.no_next_external()) {
            navigation.abort("Cannot initialize state: No external surfaces");
            return;
        }

        navigation.update_candidate(!navigation.is_init());
        update(track, navigation, cfg, ctx);

        DETRAY_VERBOSE_HOST_DEVICE("Init complete");
        DETRAY_DEBUG_HOST(detray::navigation::print_state(navigation)
                          << detray::navigation::print_candidates(
                                 navigation, cfg, track.pos(), track.dir()));
    }

    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool update(
        const track_t &track, state &navigation, const navigation::config &cfg,
        const context_type &ctx = {},
        const bool is_before_actor_run = true) const {

        DETRAY_VERBOSE_HOST_DEVICE("Called 'update()'");

        assert(navigation.is_alive());
        assert(!track.is_invalid());

        if (navigation.no_next_external()) {
            DETRAY_VERBOSE_HOST_DEVICE("No next external: Exit");
            navigation.exit();
            return false;
        }

        assert(!navigation.get_target_barcode().is_invalid());
        update_intersection(track, navigation, cfg, ctx);

        if (is_before_actor_run) {
            if (navigation.has_reached_candidate(navigation.target(), cfg)) {
                navigation.status((navigation.target().sf_desc.is_portal())
                                      ? navigation::status::e_on_portal
                                      : navigation::status::e_on_object);
                navigation.next_external();
                navigation.update_candidate(true);
                assert(navigation.has_reached_candidate(navigation.current(),
                                                        cfg));

                if (!navigation.no_next_external()) {
                    update_intersection(track, navigation, cfg, ctx);
                }

                DETRAY_VERBOSE_HOST_DEVICE("Update complete: On surface");

                // Return true to reset the step size
                return true;
            }
        }

        // Otherwise the track is moving towards a surface
        navigation.status(navigation::status::e_towards_object);

        DETRAY_VERBOSE_HOST_DEVICE("Update complete: towards surface");
        DETRAY_DEBUG_HOST(detray::navigation::print_state(navigation)
                          << detray::navigation::print_candidates(
                                 navigation, cfg, track.pos(), track.dir()));

        // Return false to scale the step size with RK4
        return false;
    }

    private:
    template <typename track_t>
    DETRAY_HOST_DEVICE inline void update_intersection(
        const track_t &track, state &navigation, const navigation::config &cfg,
        const context_type &ctx = {}) const {

        if (navigation.target().sf_desc.barcode().is_invalid()) {
            return;
        }

        const auto &det = navigation.detector();
        const auto sf = tracking_surface{det, navigation.target().sf_desc};

        const bool res =
            sf.template visit_mask<intersection_update<ray_intersector>>(
                detail::ray<algebra_type>(
                    track.pos(),
                    static_cast<scalar_type>(navigation.direction()) *
                        track.dir()),
                navigation.target(), det.transform_store(), ctx,
                cfg.template mask_tolerance<scalar_type>(),
                static_cast<scalar_type>(cfg.mask_tolerance_scalor),
                scalar_type{0.f},
                static_cast<scalar_type>(cfg.overstep_tolerance));

        // If an intersection is not found, proceed the track with safe step
        // size
        if (!res) {
            const auto path = navigation.target().path();
            navigation.update_candidate(false);
            navigation.target().set_path(
                math::copysign(navigation.safe_step_size(), path));
        }
    }
};

}  // namespace detray
