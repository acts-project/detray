/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/detector_metadata.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/benchmark/benchmark_base.hpp"
#include "tests/common/test_base/toy_detector_fixture.hpp"

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

/// Precalculate tracks
template <typename track_generator_t, typename transform3_t>
inline void fill_tracks(
    std::vector<free_track_parameters<transform3_t>> &tracks,
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
          template <typename, typename> class navigator_t = navigator,
          typename det_fixture_t = toy_detector_fixture>
struct rkn_propagation_bm : public benchmark_base {
    /// Detector dependent types
    using detector_t = detector<detector_registry::toy_detector>;
    using transform3_t = typename detector_t::transform3;
    using bfield_t = typename detector_t::bfield_type;
    using stepper_t = rk_stepper<typename bfield_t::view_t, transform3_t,
                                 unconstrained_step, stepper_policy_t>;

    /// Prefix for the benchmark name
    inline static const std::string s_name{"RKN_propagation"};

    /// Local configuration type
    struct configuration : public detray::benchmark_base::configuration {
        /// Fixture that holds the detector data
        det_fixture_t m_detector_fixture;

        /// Default construciton
        configuration() = default;

        /// Construct from a base configuration
        configuration(const detray::benchmark_base::configuration &cfg)
            : detray::benchmark_base::configuration(cfg) {}

        /// Getters
        det_fixture_t &detector_fixture() { return m_detector_fixture; }
        const det_fixture_t &detector_fixture() const {
            return m_detector_fixture;
        }
    };

    /// The benchmark configuration
    configuration m_cfg{};

    /// Default construction
    rkn_propagation_bm() = default;

    /// Construct from an externally provided configuration @param cfg
    rkn_propagation_bm(detray::benchmark_base::configuration cfg)
        : m_cfg{cfg} {}

    /// Construct from an externally provided configuration @param cfg
    rkn_propagation_bm(configuration cfg) : m_cfg{cfg} {}

    /// @return the benchmark configuration
    configuration &config() { return m_cfg; }

    std::string name() const override { return rkn_propagation_bm<>::s_name; }

    /// Prepare data and run benchmark loop
    void operator()(
        ::benchmark::State &state,
        const std::vector<free_track_parameters<transform3_t>> &tracks) const {

        using propagator_t =
            propagator<stepper_t,
                       navigator_t<detector_t, navigation::void_inspector>,
                       actor_chain<>>;

        const std::size_t n_samples{this->m_cfg.n_samples()};
        const std::size_t n_warmup{this->m_cfg.n_warmup()};

        assert(n_samples + n_warmup <= tracks.size());

        // Create detector
        vecmem::host_memory_resource host_mr;
        detector_t det =
            this->m_cfg.detector_fixture().build_detector(&host_mr);

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
                    ::benchmark::DoNotOptimize(p.propagate(p_state));
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