/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/benchmark/benchmark_base.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <random>
#include <string>
#include <thread>

namespace detray {

/// @note __plugin has to be defined with a preprocessor command
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;
using transform3 = __plugin::transform3<detray::scalar>;

/// Precalculate tracks
template <typename track_generator_t>
inline void fill_tracks(std::vector<free_track_parameters<transform3>> &tracks,
                        const typename track_generator_t::configuration &cfg) {

    // Iterate through uniformly distributed momentum directions
    track_generator_t trk_gen{cfg};
    tracks.reserve(trk_gen.size());
    for (auto traj : trk_gen) {
        tracks.push_back(traj);
    }

    // random ordering
    auto rng = std::default_random_engine{};
    std::shuffle(std::begin(tracks), std::end(tracks), rng);
}

/*template <typename bfield_bknd_t>*/
template <typename stepper_policy_t = stepper_default_policy,
          template <typename> class navigator_t = navigator>
struct rkn_toy_bm : public benchmark_base {
    /// Detector dependent types
    using detector_t = detector<detector_registry::toy_detector>;
    using transform3_t = typename detector_t::transform3;
    using bfield_t = typename detector_t::bfield_type;
    using stepper_t = rk_stepper<typename bfield_t::view_t, transform3_t,
                                 unconstrained_step, stepper_policy_t>;

    /// Prefix for the benchmark name
    inline static const std::string s_name{"RKN_toygeo"};

    /// Local configuration type
    struct configuration : public detray::benchmark_base::configuration {
        /// Number of pixel barrel layers
        unsigned int m_n_brl_layers{4u};
        /// Number of pixel endcap layers on either side
        unsigned int m_n_edc_layers{3u};

        /// Default construciton
        configuration() = default;

        /// Construct from a base configuration
        configuration(const detray::benchmark_base::configuration &cfg)
            : detray::benchmark_base::configuration(cfg) {}

        configuration &n_barrel_layers(unsigned int n) {
            m_n_brl_layers = n;
            return *this;
        }

        configuration &n_endcap_layers(unsigned int n) {
            m_n_edc_layers = n;
            return *this;
        }

        unsigned int n_barrel_layers() const { return m_n_brl_layers; }
        unsigned int n_endcap_layers() const { return m_n_edc_layers; }
    };

    /// The benchmark configuration
    configuration m_cfg{};

    /// Default construction
    rkn_toy_bm() = default;

    /// Construct from an externally provided configuration @param cfg
    rkn_toy_bm(detray::benchmark_base::configuration cfg) : m_cfg{cfg} {}

    /// Construct from an externally provided configuration @param cfg
    rkn_toy_bm(configuration cfg) : m_cfg{cfg} {}

    /// @return the benchmark configuration
    configuration &config() { return m_cfg; }

    std::string name() const override { return rkn_toy_bm<>::s_name; }

    /// Prepare data and run benchmark loop
    void operator()(
        ::benchmark::State &state,
        const std::vector<free_track_parameters<transform3_t>> &tracks) const {

        using propagator_t =
            propagator<stepper_t, navigator_t<detector_t>, actor_chain<>>;

        const std::size_t n_samples{this->m_cfg.n_samples()};
        const std::size_t n_warmup{this->m_cfg.n_warmup()};

        assert(n_samples + n_warmup <= tracks.size());

        // Create detector
        vecmem::host_memory_resource host_mr;
        detector_t det = create_toy_geometry(host_mr, m_cfg.n_barrel_layers(),
                                             m_cfg.n_endcap_layers());

        // Create propagator - should be cheap
        propagator_t p({}, {});

        // Spin down after data preparation
        if (m_cfg.do_sleep()) {
            std::this_thread::sleep_for(std::chrono::seconds(10));
        }

        // Run the benchmark
        for (auto _ : state) {
            // Warm-up
            state.PauseTiming();
            if (m_cfg.do_warmup()) {
#pragma omp parallel for
                for (unsigned int i = 0u; i < n_warmup; ++i) {
                    typename propagator_t::state p_state(tracks[i],
                                                         det.get_bfield(), det);
                    p.propagate(p_state);
                }
            }
            state.ResumeTiming();
#pragma omp parallel for
            for (unsigned int i = n_warmup; i < n_samples + n_warmup; ++i) {
                typename propagator_t::state p_state(tracks[i],
                                                     det.get_bfield(), det);
                ::benchmark::DoNotOptimize(p.propagate(p_state));
            }
        }
    }
};

}  // namespace detray