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
#include "detray/tracks/ray.hpp"
#include "detray/utils/ranges.hpp"

namespace detray {

namespace navigation {

static constexpr std::size_t default_sequence_size{100u};

}

template <typename detector_t,
          typename inspector_t = navigation::void_inspector,
          typename intersection_t =
              intersection2D<typename detector_t::surface_type,
                             typename detector_t::algebra_type, false>>
class direct_navigator {

    public:
    class state : public detray::ranges::view_interface<state> {

        // Allow the filling/updating of candidates
        friend struct intersection_initialize<ray_intersector>;
        friend struct intersection_update<ray_intersector>;

        using sequence_t = vecmem::device_vector<detray::geometry::barcode>;
        using candidate_t = intersection_type;

        public:
        using detector_type = navigator::detector_type;

        state() = delete;

        template <typename vector_t>
        DETRAY_HOST_DEVICE explicit state(const detector_type &det,
                                          const vector_t &sequence)
            : m_detector(&det) {

            // TODO: Copy sequence
            m_it = m_sequence.begin();
            m_it_rev = m_sequence.rbegin();
        }

        /// @returns the navigation heartbeat
        DETRAY_HOST_DEVICE
        bool is_alive() const { return m_heartbeat; }

        DETRAY_HOST_DEVICE
        inline auto candidate() -> candidate_t & {
            return m_candidate;
        }

        /// @return the reference of track state pointed by the iterator
        TRACCC_HOST_DEVICE
        const detray::geometry::barcode& operator()() const {
            if (m_direction == navigation::direction::e_forward) {
                return *m_it;
            } else {
                return *m_it_rev;
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
                m_it == m_track_states.end()) {
                return true;
            } else if ((m_direction != navigation::direction::e_forward) &&
                       m_it_rev == m_track_states.rend()) {
                return true;
            }
            return false;
        }

        /// @returns current navigation direction - const
        DETRAY_HOST_DEVICE
        inline auto direction() const -> navigation::direction {
            return m_direction;
        }

        /// Intersection candidate
        /// @TODO: Make it array with multiple candidates, and resort the target
        /// elements based on the path length
        candidate_t m_candidate;

        /// Target surfaces
        sequenct_t m_sequence;

        // iterator for forward direction
        typename sequence_t::iterator m_it;

        // iterator for backward direction
        typename sequence_t::reverse_iterator m_it_rev;

        /// The navigation direction
        navigation::direction m_direction{navigation::direction::e_forward};

        /// Heartbeat of this navigation flow signals navigation is alive
        bool m_heartbeat{false};
    };

    template <typename track_t>
    DETRAY_HOST_DEVICE inline void init(const track_t &track, state &navigation,
                                        const navigation::config &cfg,
                                        const context_type &ctx) const {
        navigation.m_heartbeat = true;

        auto& cand = navigation.candidate();
        cand.sf_desc.set_barcode(navigation());

        return;
    }

    template <typename track_t>
    DETRAY_HOST_DEVICE inline bool update(const track_t &track,
                                          state &navigation,
                                          const navigation::config &cfg,
                                          const context_type &ctx = {}) const {

        if (navigation.candidate().is_invalid()) {
            return false;
        }

        /*
        // Check whether this candidate is reachable by the track
        return sf.template visit_mask<intersection_update<ray_intersector>>(
            detail::ray<algebra_type>(
                track.pos(), static_cast<scalar_type>(nav_dir) * track.dir()),
            candidate, det.transform_store(), ctx,
            sf.is_portal() ? darray<scalar_type, 2>{0.f, 0.f}
                           : darray<scalar_type, 2>{cfg.min_mask_tolerance,
                                                    cfg.max_mask_tolerance},
            static_cast<scalar_type>(cfg.mask_tolerance_scalor),
            static_cast<scalar_type>(cfg.overstep_tolerance));
        */
    }
};

}  // namespace detray