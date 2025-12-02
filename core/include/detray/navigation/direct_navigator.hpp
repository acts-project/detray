/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/navigation.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/detail/intersection_kernel.hpp"
#include "detray/navigation/detail/navigation_functions.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/intersection_config.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/navigation/navigation_state.hpp"

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
        friend struct detail::intersection_update<ray_intersector>;

        template <typename state_t>
        friend constexpr void navigation::update_status(
            state_t &, const navigation::config &);

        using base_type = navigation::base_state<state, detector_type, 2u,
                                                 navigation::void_inspector,
                                                 intersection_type>;

        using candidate_t = intersection_type;
        using sf_decriptor_t = typename detector_type::surface_type;
        using nav_link_t = typename sf_decriptor_t::navigation_link;

        public:
        using value_type = candidate_t;
        using sequence_type = vecmem::device_vector<sf_decriptor_t>;

        using view_type = dvector_view<sf_decriptor_t>;
        using const_view_type = dvector_view<const sf_decriptor_t>;

        state() = delete;

        DETRAY_HOST_DEVICE constexpr state(const detector_t &det,
                                           const view_type &sequence)
            : base_type(det), m_sequence(sequence) {

            m_it = m_sequence.cbegin();
            m_it_rev = m_sequence.crbegin();

            // The next candidate is always stored in the second cache entry
            this->next_index(1u);
            this->last_index(1u);

            // Update the target with the next external surface
            set_external_target();

            assert(!m_sequence.empty());
            assert(has_next_external());
        }

        /// @returns the direct navigator always has only one candidate
        DETRAY_HOST_DEVICE
        constexpr auto n_candidates() const -> dindex { return 1u; }

        /// @returns the externally provided mask tolerance - const
        DETRAY_HOST_DEVICE
        constexpr scalar_type external_tol() const {
            return std::numeric_limits<scalar_type>::max();
        }

        /// Advance the iterator
        DETRAY_HOST_DEVICE
        constexpr sf_decriptor_t next_external() {
            assert(has_next_external());
            return is_forward() ? *m_it : *m_it_rev;
        }

        /// Set direction in which the navigator should search for candidates
        DETRAY_HOST_DEVICE
        constexpr void set_direction(const navigation::direction dir) {
            base_type::set_direction(dir);
            // Get the correct external surface from the beginning of end of
            // the sequence
            set_external_target();
        }

        /// Advance the iterator (navigation status needs to be correct)
        DETRAY_HOST_DEVICE
        constexpr void advance() {
            // The target has become the current candidate
            this->candidates()[0] = this->target();

            assert(has_next_external());
            if (is_forward()) {
                m_it++;
            } else {
                m_it_rev++;
            }

            // Update the target with the next external surface
            set_external_target();

            assert(this->target().sf_desc.is_sensitive() ||
                   this->target().sf_desc.has_material());
        }

        /// @return true if the iterator reaches the end of vector
        DETRAY_HOST_DEVICE
        constexpr bool has_next_external() const {
            return (is_forward() && m_it != m_sequence.cend()) ||
                   (!is_forward() && m_it_rev != m_sequence.crend());
        }

        /// Clear the state
        DETRAY_HOST_DEVICE constexpr void clear_cache() {
            base_type::clear_cache();
            this->next_index(1);
            this->last_index(1);
        }

        /// @returns flag that indicates whether navigation was successful
        DETRAY_HOST_DEVICE
        constexpr bool finished() const {
            // Normal exit for this navigation?
            bool is_finished = base_type::finished();

            // All external surfaces were visited?
            is_finished &= ((is_forward() && m_it == m_sequence.cend()) ||
                            (!is_forward() && m_it_rev == m_sequence.crend()));

            return is_finished;
        }

        private:
        /// @return 'true' if the navigation direction is 'forward'
        DETRAY_HOST_DEVICE
        constexpr bool is_forward() const {
            return this->direction() == navigation::direction::e_forward;
        }

        /// Set the current next external surface, depending on whether the
        /// direction is backward or forward
        DETRAY_HOST_DEVICE
        constexpr void set_external_target() {
            this->target().sf_desc = is_forward() ? *m_it : *m_it_rev;
        }

        /// Target surfaces
        sequence_type m_sequence;

        // iterator for forward direction
        typename sequence_type::const_iterator m_it;

        // iterator for backward direction
        typename sequence_type::const_reverse_iterator m_it_rev;
    };

    template <typename track_t>
    DETRAY_HOST_DEVICE inline void init(const track_t &track, state &navigation,
                                        const navigation::config &cfg,
                                        const context_type &ctx) const {
        DETRAY_VERBOSE_HOST_DEVICE("Called 'init()':");

        assert(navigation.has_next_external());

        // Clean up state
        navigation.clear_cache();

        // Update the next candidate
        update(track, navigation, cfg, ctx);

        DETRAY_VERBOSE_HOST_DEVICE("Init complete!");
    }

    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool update(
        const track_t &track, state &navigation, const navigation::config &cfg,
        const context_type &ctx = {},
        const bool /*is_before_actor_run*/ = true) const {

        DETRAY_VERBOSE_HOST_DEVICE("Called 'update()'");

        // Do not resurrect a failed/finished navigation state
        assert(navigation.is_alive());
        assert(!track.is_invalid());

        if (!navigation.has_next_external()) {
            DETRAY_VERBOSE_HOST_DEVICE("No next external: Exit");
            navigation.exit();
            return false;
        }

        const detector_type &det = navigation.detector();
        constexpr bool is_init{true};

        // Update only the current candidate and the corresponding next target
        if (navigation.trust_level() != navigation::trust_level::e_full) {
            if (!navigation::update_candidate(
                    navigation.direction(), navigation.target(), track, det,
                    cfg.intersection, navigation.external_tol(), ctx)) {

                DETRAY_DEBUG_HOST("Update candidate: After "
                                  << navigation.target());
                navigation.abort("Could not find next surface");
                return false;
            }

            // Update navigation flow on the new candidate information and set
            // the next target if the surface has been reached
            navigation::update_status(navigation, cfg);

            if (navigation.trust_level() != navigation::trust_level::e_full) {
                DETRAY_DEBUG_HOST("OH NO!");
            }

            assert(navigation.trust_level() == navigation::trust_level::e_full);

            // The work is done if: the track has not reached a surface yet
            // or trust is gone (portal was reached or the cache is broken).
            if (navigation.status() == navigation::status::e_towards_object) {

                DETRAY_VERBOSE_HOST_DEVICE("Update complete (towards object):");
                DETRAY_DEBUG_HOST("\n"
                                  << detray::navigation::print_state(navigation)
                                  << detray::navigation::print_candidates(
                                         navigation, cfg, track.pos(),
                                         track.dir()));

                return !is_init;
            }

            // Track is on surface: Update the new target (if available)
            if (navigation.has_next_external() &&
                !navigation::update_candidate(
                    navigation.direction(), navigation.target(), track, det,
                    cfg.intersection, navigation.external_tol(), ctx)) {
                DETRAY_DEBUG_HOST("Update candidate: After "
                                  << navigation.target());

                navigation.abort("Could not find next surface");
                return false;
            }

            // At this point, the track has to be on surface:
            // Set volume index to the next volume provided by the portal
            navigation.set_volume(navigation.current().sf_desc.volume());
            navigation.trust_level(navigation::trust_level::e_full);

            DETRAY_VERBOSE_HOST_DEVICE("Update complete (on surface):");
            DETRAY_DEBUG_HOST("\n"
                              << detray::navigation::print_state(navigation)
                              << detray::navigation::print_candidates(
                                     navigation, cfg, track.pos(),
                                     track.dir()));

            // Return true to reset the step size of the RKN algorithm
            return is_init;
        }

        DETRAY_VERBOSE_HOST_DEVICE(" -> Full trust: Nothing left to do");

        return !is_init;
    }
};

}  // namespace detray
