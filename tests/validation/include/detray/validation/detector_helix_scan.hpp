/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/types.hpp"
#include "tests/common/test_base/fixture_base.hpp"
#include "tests/common/tools/particle_gun.hpp"
#include "tests/common/tools/ray_scan_utils.hpp"

// System include(s)
#include <iostream>
#include <string>

namespace detray {

/// @brief Test class that runs the helix scan on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class helix_scan : public test::fixture_base<> {

    using transform3_t = typename detector_t::transform3;
    using free_track_parameters_t = free_track_parameters<transform3_t>;

    public:
    using fixture_type = test::fixture_base<>;

    struct config : public fixture_type::configuration {
        using trk_gen_config_t = typename uniform_track_generator<
            free_track_parameters_t>::configuration;

        std::string m_name{"helix_scan"};
        bool m_write_inters{false};
        trk_gen_config_t m_trk_gen_cfg{};

        /// Getters
        /// @{
        const std::string &name() const { return m_name; }
        bool write_intersections() const { return m_write_inters; }
        trk_gen_config_t &track_generator() { return m_trk_gen_cfg; }
        const trk_gen_config_t &track_generator() const {
            return m_trk_gen_cfg;
        }
        /// @}

        /// Setters
        /// @{
        config &name(const std::string n) {
            m_name = n;
            return *this;
        }
        config &write_intersections(const bool do_write) {
            m_write_inters = do_write;
            return *this;
        }
        /// @}
    };

    template <typename config_t>
    explicit helix_scan(const detector_t &det,
                        const typename detector_t::name_map &names,
                        const config_t &cfg = {})
        : m_det{det}, m_names{names} {
        m_cfg.name(cfg.name());
        m_cfg.write_intersections(cfg.write_intersections());
        m_cfg.overstepping_tolerance(cfg.overstepping_tolerance());
        m_cfg.track_generator() = cfg.track_generator();
    }

    /// Run the helix scan
    void TestBody() override {
        using b_field_t = typename detector_t::bfield_type;
        using scalar_t = typename detector_t::scalar_type;
        using nav_link_t = typename detector_t::surface_type::navigation_link;

        constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};

        // Index of the volume that the helix origin lies in
        dindex start_index{0u};

        // B-field vector for helix
        const auto &origin = m_cfg.track_generator().origin();
        typename b_field_t::view_t bfield_view{m_det.get_bfield()};
        const typename b_field_t::output_t bvec =
            bfield_view.at(origin[0], origin[1], origin[2]);
        const typename fixture_type::point3 B{bvec[0], bvec[1], bvec[2]};

        // Iterate through uniformly distributed momentum directions
        std::size_t n_tracks{0u};
        auto trk_state_generator =
            uniform_track_generator<free_track_parameters_t>(
                m_cfg.track_generator());

        std::cout << "INFO: Running helix scan on: " << m_names.at(0) << "\n("
                  << trk_state_generator.size() << " helices) ...\n"
                  << std::endl;

        for (auto trk : trk_state_generator) {

            // Prepare for overstepping in the presence of b fields
            trk.set_overstep_tolerance(m_cfg.overstepping_tolerance());

            // Get ground truth helix from track
            detail::helix helix(trk, &B);

            // Shoot helix through the detector and record all surfaces it
            // encounters
            const auto intersection_record = particle_gun::shoot_particle(
                m_det, helix, 14.9f * unit<scalar_t>::um);

            // Create a trace of the volume indices that were encountered
            // and check that portal intersections are connected
            auto [portal_trace, surface_trace] =
                trace_intersections<leaving_world>(intersection_record,
                                                   start_index);

            // Is the succession of volumes consistent ?
            ASSERT_TRUE(check_connectivity<leaving_world>(portal_trace))
                << "\nFailed on helix " << n_tracks << "/"
                << trk_state_generator.size() << "\n"
                << helix;

            ++n_tracks;
        }
    }

    private:
    /// The configuration of this test
    config m_cfg;
    /// The detector to be checked
    const detector_t &m_det;
    /// Volume names
    const typename detector_t::name_map &m_names;
};

}  // namespace detray
