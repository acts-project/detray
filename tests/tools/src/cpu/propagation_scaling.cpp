/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/detectors/bfield.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/type_list.hpp"

// Detray IO include(s)
#include "detray/io/frontend/detector_reader.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/cpu/propagation_benchmark.hpp"

// Detray test include(s).
#include "detray/test/utils/simulation/event_generator/track_generators.hpp"
#include "detray/test/utils/types.hpp"

// Detray tools include(s)
#include "detray/options/detector_io_options.hpp"
#include "detray/options/parse_options.hpp"
#include "detray/options/propagation_options.hpp"
#include "detray/options/track_generator_options.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <algorithm>
#include <string>
#include <thread>

namespace po = boost::program_options;

using namespace detray;

int main(int argc, char** argv) {

    // Use the most general type to be able to read in all detector files
    using detector_t = detray::detector<test::default_metadata>;
    using test_algebra = typename detector_t::algebra_type;
    using scalar = dscalar<test_algebra>;
    using vector3 = dvector3D<test_algebra>;

    using free_track_parameters_t = free_track_parameters<test_algebra>;
    using uniform_gen_t =
        detail::random_numbers<scalar, std::uniform_real_distribution<scalar>>;
    using track_generator_t =
        random_track_generator<free_track_parameters_t, uniform_gen_t>;

    using field_t = bfield::const_field_t<scalar>;
    using stepper_t = rk_stepper<typename field_t::view_t, test_algebra>;
    using empty_chain_t = actor_chain<>;
    using default_chain =
        actor_chain<parameter_transporter<test_algebra>,
                    pointwise_material_interactor<test_algebra>,
                    parameter_resetter<test_algebra>>;

    // Code for the openMP scheduling scheme,
    // @see https://www.openmp.org/spec-html/5.0/openmpsu121.html
    // omp_sched_static = 0x1,
    // omp_sched_dynamic = 0x2,
    // omp_sched_guided = 0x3,
    // omp_sched_auto = 0x4,
    //
    // schedule modifier
    // omp_sched_monotonic = 0x80000000u
    constexpr int sched_policy{1};
    /// Number of host threads
    std::vector<int> n_threads{1, 2, 4, 8, 16, 32, 64, 128, 256};

    /// Number of tracks per thread in weak scaling
    constexpr int chunk_size{1000};
    std::vector<int> n_tracks_weak_sc;
    // Ensure that number of tracks is divisible by number of threads
    for (const int n : n_threads) {
        n_tracks_weak_sc.push_back(n * chunk_size);
    }

    /// Maximum number of tracks to be assigned to thread (invlaid:
    /// set per benchmark case as strong_sc_sample_size/#threads)
    constexpr int max_chunk_size{detail::invalid_value<int>()};
    /// Strong scaling sample size
    constexpr std::size_t strong_sc_sample_size{51'200u};

    // Host memory resource
    vecmem::host_memory_resource host_mr;

    // Constant magnetic field
    vector3 B{0.f, 0.f, 2.f * unit<scalar>::T};

    //
    // Configuration
    //

    // Google benchmark specific options
    ::benchmark::Initialize(&argc, argv);

    // Specific options for this test
    po::options_description desc("\ndetray propagation scaling options");

    desc.add_options()("context", po::value<dindex>(),
                       "Index of the geometry context")(
        "bknd_name", po::value<std::string>(), "Name of the Processor")(
        "sort_tracks", "Does nothing in scaling test");

    // Configs to be filled
    detray::io::detector_reader_config reader_cfg{};
    track_generator_t::configuration trk_cfg{};
    propagation::config prop_cfg{};
    detray::benchmarks::benchmark_base::configuration bench_cfg{};

    // Read options from commandline
    po::variables_map vm = detray::options::parse_options(
        desc, argc, argv, reader_cfg, trk_cfg, prop_cfg);

    // The geometry context to be used
    detector_t::geometry_context gctx;
    if (vm.count("context")) {
        gctx = detector_t::geometry_context{vm["context"].as<dindex>()};
    }
    std::string proc_name{"unknown"};
    if (vm.count("bknd_name")) {
        proc_name = vm["bknd_name"].as<std::string>();
    }

    //
    // Prepare data
    //

    // Read the detector geometry
    reader_cfg.do_check(true);

    const auto [det, names] =
        detray::io::read_detector<detector_t>(host_mr, reader_cfg);
    const std::string& det_name = det.name(names);

    // Generate the track samples for weak scaling
    auto track_samples_weak_sc =
        detray::benchmarks::generate_track_samples<track_generator_t>(
            &host_mr, n_tracks_weak_sc, trk_cfg, false);

    // Generate track sample for strong scaling
    trk_cfg.n_tracks(strong_sc_sample_size);
    auto single_sample = detray::benchmarks::generate_tracks<track_generator_t>(
        &host_mr, trk_cfg, false);
    std::vector<dvector<free_track_parameters_t>> track_samples_strong_sc{
        std::move(single_sample)};

    // Create a constant b-field
    auto bfield = bfield::create_const_field<scalar>(B);

    // Build actor states
    dtuple<> empty_state{};

    parameter_transporter<test_algebra>::state transporter_state{};
    pointwise_material_interactor<test_algebra>::state interactor_state{};
    parameter_resetter<test_algebra>::state resetter_state{};

    auto actor_states = detail::make_tuple<dtuple>(
        transporter_state, interactor_state, resetter_state);

    //
    // Register benchmarks
    //

    // Number of warmup tracks
    const int n_max_tracks{*std::ranges::max_element(n_tracks_weak_sc)};
    bench_cfg.n_warmup(
        static_cast<int>(std::ceil(0.1f * static_cast<float>(n_max_tracks))));

    if (prop_cfg.stepping.do_covariance_transport) {
        // Number of tracks to be sampled and number of threads are the same
        detray::benchmarks::register_benchmark<
            detray::benchmarks::host_propagation_bm, stepper_t, default_chain>(
            det_name + "_W_COV_TRANSPORT_WEAK-SCALING", bench_cfg, prop_cfg,
            det, bfield, &actor_states, track_samples_weak_sc, n_tracks_weak_sc,
            n_threads, chunk_size, sched_policy);

        // Number of tracks is increased on a fixed size sample
        detray::benchmarks::register_benchmark<
            detray::benchmarks::host_propagation_bm, stepper_t, default_chain>(
            det_name + "_W_COV_TRANSPORT_STRONG-SCALING", bench_cfg, prop_cfg,
            det, bfield, &actor_states, track_samples_strong_sc,
            {strong_sc_sample_size}, n_threads, max_chunk_size, sched_policy);
    } else {
        detray::benchmarks::register_benchmark<
            detray::benchmarks::host_propagation_bm, stepper_t, empty_chain_t>(
            det_name + "_WEAK-SCALING", bench_cfg, prop_cfg, det, bfield,
            &empty_state, track_samples_weak_sc, n_tracks_weak_sc, n_threads,
            chunk_size, sched_policy);

        detray::benchmarks::register_benchmark<
            detray::benchmarks::host_propagation_bm, stepper_t, empty_chain_t>(
            det_name + "_STRONG-SCALING", bench_cfg, prop_cfg, det, bfield,
            &empty_state, track_samples_strong_sc, {strong_sc_sample_size},
            n_threads, max_chunk_size, sched_policy);
    }

    // These fields are needed by the plotting scripts
    ::benchmark::AddCustomContext("Backend", "CPU");
    ::benchmark::AddCustomContext("Backend Name", proc_name);
    ::benchmark::AddCustomContext(
        "Max no. Threads", std::to_string(std::thread::hardware_concurrency()));
    ::benchmark::AddCustomContext("Algebra-plugin",
                                  detray::types::get_name<test_algebra>());

    // Run benchmarks
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
}
