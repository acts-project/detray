/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/benchmarks/device/cuda/propagator.hpp"
#include "detray/definitions/detail/cuda_definitions.hpp"
#include "detray/utils/logging.hpp"

// CUDA include(s)
#include <cooperative_groups.h>

#include <cuda/pipeline>

// Disables `pipeline_shared_state` initialization warning.
#pragma nv_diag_suppress static_var_with_dynamic_init

namespace detray {

/// Run the navigator
template <typename navigator_t>
__device__ inline void navigate(
    const navigation::config &cfg,
    const typename navigator_t::detector_type::geometry_context &gtx,
    stepping::result<typename navigator_t::detector_type::algebra_type, 32u>
        &stepper_result,
    typename navigator_t::state &navigation) {

    using detector_t = typename navigator_t::detector_type;

    auto block = cooperative_groups::this_thread_block();
    const unsigned int wid{block.thread_rank() % 32u};

    assert(wid < 32u);

    constexpr navigator_t navigator{};

    navigator.update(stepper_result.free_param(wid), navigation, cfg, gtx);
}

/// Run the stepper
template <typename stepper_t, typename actor_chain_t, typename metadata_t>
__device__ inline void take_step(
    const stepping::config &cfg,
    navigation::result<metadata_t, 32u> &nav_result,
    typename stepper_t::state &step_state,
    typename actor_chain_t::state_ref_tuple &actor_states) {

    using algebra_t = typename metadata_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;

    auto block = cooperative_groups::this_thread_block();
    const unsigned int wid{block.thread_rank() % 32u};

    assert(wid < 32u);

    constexpr stepper_t stepper{};

    // Set access to the volume material for the stepper
    /*auto vol = nav_result.current_volume(wid);
    const material<scalar_t> *vol_mat_ptr =
        vol.has_material() ? vol.material_parameters(step_state().pos())
                           : nullptr;*/

    // Break automatic step size scaling by the stepper when a
    // surface was reached and whenever the navigation is
    // (re-)initialized
    const bool reset_stepsize{nav_result.is_on_surface(wid) || true};

    // Take the step
    stepper.step(nav_result(wid), step_state, cfg, reset_stepsize, nullptr);

    // Reduce navigation trust level according to stepper update
    // typename stepper_t::policy_type{}(stepping.policy_state(),
    //                                    propagation);
    nav_result.set_high_trust(wid);
}

/// Run the actors
template <typename stepper_t, typename actor_chain_t, typename metadata_t>
__device__ inline void run_actors(
    navigation::result<metadata_t, 32u> &nav_result,
    typename stepper_t::state &step_state,
    typename actor_chain_t::state_ref_tuple &actor_states) {

    using stepper_state_t = typename stepper_t::state;
    using actor_states_t = typename actor_chain_t::state_ref_tuple;

    constexpr actor_chain_t run_actors{};
}

/// Specialize the warps to run stepping, navigation and actors independently
template <typename navigator_t, typename stepper_t, typename actor_chain_t>
__global__ void __launch_bounds__(256, 4) propagation_kernel(
    const __grid_constant__ propagation::config cfg,
    const __grid_constant__
    typename navigator_t::detector_type *const pinned_detector_ptr,
    const __grid_constant__ typename stepper_t::magnetic_field_type field_view,
    vecmem::data::vector_view<typename stepper_t::free_track_parameters_type>
        tracks_view,
    vecmem::data::vector_view<typename navigator_t::state>
        navigator_states_view,
    vecmem::data::vector_view<typename stepper_t::state> stepper_states_view,
    const __grid_constant__
    typename actor_chain_t::state_tuple *const pinned_actor_state_ptr) {

    __syncthreads();
    cuda::std::chrono::high_resolution_clock::time_point block_t1 =
        cuda::std::chrono::high_resolution_clock::now();

    using detector_device_t = typename navigator_t::detector_type;
    using metadata_t = typename detector_device_t::metadata;
    using algebra_t = typename detector_device_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using nav_link_t =
        typename detector_device_t::surface_type::navigation_link;

    using track_t = typename stepper_t::free_track_parameters_type;
    using stepper_state_t = typename stepper_t::state;
    using navigation_state_t = typename navigator_t::state;
    using stepper_result_t = stepping::result<algebra_t, 32u>;
    using navigation_result_t = navigation::result<metadata_t, 32u>;

    int nav_cycles{0};
    int step_cycles{0};
    double nav_work{0.f};
    double step_work{0.f};
    double nav_mem{0.f};
    double step_mem{0.f};
    double nav_wait{0.f};
    double step_wait{0.f};
    double nav_total{0.f};
    double step_total{0.f};
    cuda::std::chrono::high_resolution_clock::time_point init_t1 =
        cuda::std::chrono::high_resolution_clock::now();

    constexpr int n_propagation_steps{100};
    constexpr unsigned int n_warps{256u / 32u};

    assert(blockDim.y == blockDim.z == 1u);
    assert(blockDim.x % 64u == 0u);
    assert(n_warps % 2u == 0u);

    auto block = cooperative_groups::this_thread_block();

    // Divide the block into warp pairs, where work is divided between
    // the two warps into nav <-> step/act
    cooperative_groups::thread_block_tile<64> warp_pair =
        cooperative_groups::tiled_partition<64>(block);

    // Thread indices
    const unsigned int gid = threadIdx.x + blockIdx.x * blockDim.x;
    const unsigned int bid = block.thread_rank();
    const unsigned int wid{bid % 32u};   //< Thread rank in warp
    const unsigned int wpid{bid / 64u};  //< Warp pair index per block

    // Geometry context
    const typename detector_device_t::geometry_context gctx{};

    // Input data:
    // One navigation and one stepping/acting thread work on 2 tracks (2:2)
    const vecmem::device_vector<track_t> tracks(tracks_view);

    // Setup the navigator and stepper states buffers to memoize states between
    // kernel calls
    vecmem::device_vector<navigation_state_t> navigator_states{
        navigator_states_view};
    vecmem::device_vector<stepper_state_t> stepper_states{stepper_states_view};

    assert(navigator_states.size() == tracks.size());
    assert(stepper_states.size() == tracks.size());

    /// The navigation and stepper result data are exchanged between warps
    __shared__ navigation_result_t navigator_results[n_warps];
    __shared__ stepper_result_t stepper_results[n_warps];

    // Create the propagation pipeline.
    constexpr auto scope{cuda::thread_scope_block};
    constexpr int n_stages{2};

    // One pipeline per warp pair: 256 / 64 = 4
    __shared__ cuda::pipeline_shared_state<scope, n_stages>
        shared_states[n_warps / 2u];

    // Print the threads indices
    DETRAY_INFO_DEVICE("Thread %d (%d) in block %d has warp pair index: %d",
                       gid, bid, blockIdx.x, wpid);

    // Split the warp pair into their roles
    const bool does_navigation{warp_pair.thread_rank() < warp_pair.size() / 2u};
    auto thread_role = does_navigation ? cuda::pipeline_role::consumer
                                       : cuda::pipeline_role::producer;

    // Create the partitioned propagation pipeline for the warp pair
    auto propagation =
        cuda::make_pipeline(warp_pair, &shared_states[wpid], thread_role);

    // Each propagation step takes four cycles through the pipelilne to complete
    constexpr int n_pipline_cycles{4 * n_propagation_steps};

    // Map the id in the worker warps back onto the same data offset
    const unsigned int data_idx{does_navigation ? gid : gid - 32u};

    // Have all threads participate in the propagation initialization:
    // The stepper warps move the result via the global navigation states to the
    // navigation warps local memory
    navigation_state_t tmp_nav_state{*pinned_detector_ptr};
    {
        // Navigation stage: which of the two tracks per warp pair to navigate
        const unsigned int stage{does_navigation ? 0u : 1u};
        // Fetch the correct track and nav/step result
        const unsigned int i{data_idx + stage * 32u};
        const unsigned int j{i / 32u - blockIdx.x * n_warps};

        if (i < tracks.size()) {
            // Read initial track data
            // stepper_states.at(i) = stepper_state_t{tracks.at(i), field_view};

            // Take the navigation state slot as temporary scratch space
            navigation_state_t &nav_state =
                does_navigation ? tmp_nav_state : navigator_states.at(i);
            nav_state = navigation_state_t{*pinned_detector_ptr};

            constexpr navigator_t navigator{};
            navigator.init(tracks.at(i), nav_state, cfg.navigation, gctx);

            if (wid == 0) {
                navigator_results[j] = {};
            }
            navigator_results[j].m_volume_index[wid] = nav_state.volume();
            navigator_results[j].m_dist_to_next[wid] = nav_state();
            navigator_results[j].m_status[wid] = nav_state.status();
            navigator_results[j].m_trust_level[wid] = nav_state.trust_level();

            DETRAY_INFO_DEVICE(
                "Thread %d: Navigation init: %f mm (stage %d, data %d, res %d)",
                block.thread_rank(), nav_state(), stage, i, j);

            assert(nav_state.trust_level() == navigation::trust_level::e_full);
            assert(nav_state.is_alive());

            if (nav_state.trust_level() != navigation::trust_level::e_full ||
                !nav_state.is_alive()) {
                DETRAY_ERROR_DEVICE(
                    "Thread %d: INIT FAILED %f mm (stage %d, data %d, res %d)",
                    block.thread_rank(), nav_state(), stage, i, j);
            }
        }
    }
    warp_pair.sync();

    cuda::std::chrono::high_resolution_clock::time_point init_t2 =
        cuda::std::chrono::high_resolution_clock::now();

    const auto init_time =
        cuda::std::chrono::duration_cast<cuda::std::chrono::duration<double>>(
            init_t2 - init_t1);

    //
    // Run propagation steps
    //

    // The first half of the warp pair performs navigation, the second half
    // stepping and actor calls
    cuda::std::chrono::high_resolution_clock::time_point prop_t1 =
        cuda::std::chrono::high_resolution_clock::now();
    if (does_navigation) {

        cuda::std::chrono::high_resolution_clock::time_point nav_t1 =
            cuda::std::chrono::high_resolution_clock::now();

        if (data_idx + 32u < tracks.size()) {
            DETRAY_INFO_DEVICE(
                "Thread %d (%d) in block %d runs navigation (tracks: %d, %d)",
                gid, bid, blockIdx.x, data_idx, data_idx + 32u);
        } else if (data_idx < tracks.size()) {
            DETRAY_INFO_DEVICE(
                "Thread %d (%d) in block %d runs navigation (tracks: %d, -)",
                gid, bid, blockIdx.x, data_idx);
        }

        // Fetch the navigation states into local memory for coalesced access
        // TODO: Remove once the memory layout of the states is fully SoA
        std::array<navigation_state_t, 2u> navigation{
            std::move(tmp_nav_state), navigator_states[data_idx + 32u]};

        if (data_idx < tracks.size()) {
            [[maybe_unused]] const auto &pos = tracks.at(data_idx).pos();
            DETRAY_INFO_DEVICE(
                "Thread %d: Nav state %f, pos [r: %f, z: %f] (stage 0, data "
                "%d)",
                block.thread_rank(), navigation[0](), vector::perp(pos), pos[2],
                data_idx);
        }
        if (data_idx + 32u < tracks.size()) {
            [[maybe_unused]] const auto &pos = tracks.at(data_idx + 32u).pos();
            DETRAY_INFO_DEVICE(
                "Thread %d: Nav state %f, pos [r: %f, z: %f] (stage 1, data "
                "%d)",
                block.thread_rank(), navigation[1](), vector::perp(pos), pos[2],
                data_idx + 32u);
        }

        // Propagation loop
        int stage{0};  //< which of the two tracks to navigate
        for (std::size_t step = 0u; step < n_pipline_cycles; ++step) {

            const unsigned int i{data_idx + stage * 32u};
            const unsigned int j{i / 32u - blockIdx.x * n_warps};

            // Wait for the stepper/actor to commit next step/action to pipeline
            //__syncwarp();
            cuda::std::chrono::high_resolution_clock::time_point wait_t1 =
                cuda::std::chrono::high_resolution_clock::now();
            propagation.consumer_wait();
            cuda::std::chrono::high_resolution_clock::time_point wait_t2 =
                cuda::std::chrono::high_resolution_clock::now();

            const auto wait_time = cuda::std::chrono::duration_cast<
                cuda::std::chrono::duration<double>>(wait_t2 - wait_t1);
            nav_wait += wait_time.count();

            cuda::std::chrono::high_resolution_clock::time_point cycle_t1 =
                cuda::std::chrono::high_resolution_clock::now();

            navigation_result_t &nav_res = navigator_results[j];
            navigation_state_t &nav_state = navigation[stage];

            if (i < tracks.size() && nav_res.is_alive(wid) &&
                nav_res.trust_level(wid) != navigation::trust_level::e_full) {
                // Make sure the 'trust level' (configured by the actors) is set
                // TODO: Add function to fully configure nav. state from result
                nav_state.set_high_trust();
                nav_cycles++;

                [[maybe_unused]] const auto &pos =
                    stepper_results[j].free_param(wid).pos();
                DETRAY_INFO_DEVICE(
                    "STEP %ld: Thread %d (Nav): pos [%f, %f, %f] [r: %f, z: "
                    "%f] (stage %d, "
                    "data %d, res %d)",
                    step - stage, block.thread_rank(), pos[0], pos[1], pos[2],
                    vector::perp(pos), pos[2], stage, i, j);

                navigate<navigator_t>(cfg.navigation, gctx, stepper_results[j],
                                      nav_state);

                if (nav_state.is_alive()) {
                    DETRAY_INFO_DEVICE(
                        "STEP %ld: Thread %d: Dist %f (stage %d, data %d, res "
                        "%d)",
                        step - stage, block.thread_rank(), nav_state(), stage,
                        i, j);
                } else if (nav_state.finished()) {
                    DETRAY_INFO_DEVICE(
                        "STEP %ld: Thread %d: Navigation finished (stage %d, "
                        "data %d, res %d)",
                        step - stage, block.thread_rank(), stage, i, j);
                } else {
                    printf(
                        "STEP %ld: Thread %d: Navigation aborted (stage %d, "
                        "data %d, res %d)",
                        step - stage, block.thread_rank(), stage, i, j);
                }

                if (nav_state.trust_level() !=
                        navigation::trust_level::e_full &&
                    nav_state.is_alive()) {
                    DETRAY_ERROR_DEVICE(
                        "STEP %ld: Thread %d: NAV FAILED %f (stage %d, data "
                        "%d, res %d)",
                        step - stage, block.thread_rank(), nav_state(), stage,
                        i, j);
                }
            }

            cuda::std::chrono::high_resolution_clock::time_point cycle_t2 =
                cuda::std::chrono::high_resolution_clock::now();

            // Publish navigator results
            //__syncwarp();

            cuda::std::chrono::high_resolution_clock::time_point mem_t1 =
                cuda::std::chrono::high_resolution_clock::now();

            // Move navigation data from global to shared memory
            cuda::memcpy_async(&(nav_res.m_volume_index[wid]),
                               &(nav_state.m_volume_index), sizeof(nav_link_t),
                               propagation);
            cuda::memcpy_async(
                &(nav_res.m_dist_to_next[wid]),
                &(nav_state
                      .m_candidates[static_cast<std::size_t>(nav_state.m_next)]
                      .ip.path),
                sizeof(scalar_t), propagation);
            cuda::memcpy_async(&(nav_res.m_status[wid]), &(nav_state.m_status),
                               sizeof(navigation::status), propagation);
            cuda::memcpy_async(&(nav_res.m_trust_level[wid]),
                               &(nav_state.m_trust_level),
                               sizeof(navigation::trust_level), propagation);

            cuda::std::chrono::high_resolution_clock::time_point mem_t2 =
                cuda::std::chrono::high_resolution_clock::now();

            const auto mem_time = cuda::std::chrono::duration_cast<
                cuda::std::chrono::duration<double>>(mem_t2 - mem_t1);

            nav_mem += mem_time.count();

            propagation.consumer_release();

            // Flip stage and navigate the other track
            stage = (stage + 1) % n_stages;

            const auto cycle_time = cuda::std::chrono::duration_cast<
                cuda::std::chrono::duration<double>>(cycle_t2 - cycle_t1);

            nav_work += cycle_time.count();
        }

        DETRAY_INFO_DEVICE("Thread %d, block %d (Nav): exit pipeline",
                           block.thread_rank(), blockIdx.x);

        cuda::std::chrono::high_resolution_clock::time_point nav_t2 =
            cuda::std::chrono::high_resolution_clock::now();

        const auto nav_time = cuda::std::chrono::duration_cast<
            cuda::std::chrono::duration<double>>(nav_t2 - nav_t1);
        nav_total += nav_time.count();

    } else {
        cuda::std::chrono::high_resolution_clock::time_point step_t1 =
            cuda::std::chrono::high_resolution_clock::now();

        if (data_idx + 32u < tracks.size()) {
            DETRAY_INFO_DEVICE(
                "Thread %d (%d) in block %d runs stepping (tracks: %d, %d)",
                gid, bid, blockIdx.x, data_idx, data_idx + 32u);
        } else if (data_idx < tracks.size()) {
            DETRAY_INFO_DEVICE(
                "Thread %d (%d) in block %d runs stepping (tracks: %d, -)", gid,
                bid, blockIdx.x, data_idx);
        }

        // Local stepper states (e.g. track position and direction, B-field)
        std::array<stepper_state_t, 2u> stepping{stepper_state_t{field_view},
                                                 stepper_state_t{field_view}};

        // Make a copy of the initial actor states (blueprint) into local
        // memory
        std::array<typename actor_chain_t::state_tuple, 2u> actor_states{
            *pinned_actor_state_ptr, *pinned_actor_state_ptr};

        using state_refs_t = typename actor_chain_t::state_ref_tuple;
        std::array<state_refs_t, 2u> actor_state_refs{
            actor_chain_t::setup_actor_states(actor_states[0]),
            actor_chain_t::setup_actor_states(actor_states[1])};

        // Initialize the propagation pipeline
        for (int stage = 0; stage < 2; ++stage) {
            const unsigned int i{data_idx + stage * 32u};
            const unsigned int j{i / 32u - blockIdx.x * n_warps};

            // __syncwarp();
            cuda::std::chrono::high_resolution_clock::time_point wait_t1 =
                cuda::std::chrono::high_resolution_clock::now();
            propagation.producer_acquire();
            cuda::std::chrono::high_resolution_clock::time_point wait_t2 =
                cuda::std::chrono::high_resolution_clock::now();

            const auto wait_time = cuda::std::chrono::duration_cast<
                cuda::std::chrono::duration<double>>(wait_t2 - wait_t1);
            step_wait += wait_time.count();

            cuda::std::chrono::high_resolution_clock::time_point cycle_t1 =
                cuda::std::chrono::high_resolution_clock::now();

            stepper_state_t &step_state = stepping[stage];

            if (i < tracks.size()) {

                step_state() = tracks.at(i);

                // Check actor update
                navigation_result_t &nav_res = navigator_results[j];
                if constexpr (std::same_as<actor_chain_t, actor_chain<>>) {
                    take_step<stepper_t, actor_chain_t>(
                        cfg.stepping, nav_res, step_state,
                        actor_state_refs[stage]);

                    DETRAY_INFO_DEVICE(
                        "Thread %d: Fill pipeline: path length %f mm "
                        "(stage "
                        "%d, data %d, "
                        "res %d)",
                        block.thread_rank(), step_state.path_length(), stage, i,
                        j);
                } else {
                    run_actors<stepper_t, actor_chain_t>(
                        nav_res, step_state, actor_state_refs[stage]);

                    DETRAY_INFO_DEVICE(
                        "Thread %d: Fill pipeline: Run actors (stage %d, "
                        "data "
                        "%d, "
                        "res %d)",
                        block.thread_rank(), stage, i, j);
                }
            }

            cuda::std::chrono::high_resolution_clock::time_point cycle_t2 =
                cuda::std::chrono::high_resolution_clock::now();

            //__syncwarp();

            cuda::std::chrono::high_resolution_clock::time_point mem_t1 =
                cuda::std::chrono::high_resolution_clock::now();

            // Publish intial stepper/actor results (global ->
            // shared)(stepper_results[j].pos0.data() + wid,
            cuda::memcpy_async(stepper_results[j].pos0.data() + wid,
                               step_state.m_track.m_vector[0].data(),
                               sizeof(scalar_t), propagation);
            cuda::memcpy_async(stepper_results[j].pos1.data() + wid,
                               step_state.m_track.m_vector[0].data() + 1,
                               sizeof(scalar_t), propagation);
            cuda::memcpy_async(stepper_results[j].pos2.data() + wid,
                               step_state.m_track.m_vector[0].data() + 2,
                               sizeof(scalar_t), propagation);
            cuda::memcpy_async(stepper_results[j].dir0.data() + wid,
                               step_state.m_track.m_vector[0].data() + 4,
                               sizeof(scalar_t), propagation);
            cuda::memcpy_async(stepper_results[j].dir1.data() + wid,
                               step_state.m_track.m_vector[0].data() + 5,
                               sizeof(scalar_t), propagation);
            cuda::memcpy_async(stepper_results[j].dir2.data() + wid,
                               step_state.m_track.m_vector[0].data() + 6,
                               sizeof(scalar_t), propagation);

            cuda::std::chrono::high_resolution_clock::time_point mem_t2 =
                cuda::std::chrono::high_resolution_clock::now();

            const auto mem_time = cuda::std::chrono::duration_cast<
                cuda::std::chrono::duration<double>>(mem_t2 - mem_t1);

            step_mem += mem_time.count();

            propagation.producer_commit();

            const auto cycle_time = cuda::std::chrono::duration_cast<
                cuda::std::chrono::duration<double>>(cycle_t2 - cycle_t1);
            step_work += cycle_time.count();
        }

        // Propagation loop
        int stage = 0;                      //< which track to advance
        std::array<int, 2u> prop_phase{0};  //< whether to run actors or stepper
        for (std::size_t step = 0; step < n_pipline_cycles; ++step) {
            const unsigned int i{data_idx + stage * 32u};
            const unsigned int j{i / 32u - blockIdx.x * n_warps};

            // Wait for navigation to finish
            cuda::std::chrono::high_resolution_clock::time_point wait_t1 =
                cuda::std::chrono::high_resolution_clock::now();
            propagation.producer_acquire();
            cuda::std::chrono::high_resolution_clock::time_point wait_t2 =
                cuda::std::chrono::high_resolution_clock::now();

            const auto wait_time = cuda::std::chrono::duration_cast<
                cuda::std::chrono::duration<double>>(wait_t2 - wait_t1);
            step_wait += wait_time.count();

            cuda::std::chrono::high_resolution_clock::time_point cycle_t1 =
                cuda::std::chrono::high_resolution_clock::now();

            navigation_result_t &nav_res = navigator_results[j];
            stepper_state_t &step_state = stepping[stage];

            if (i < tracks.size() && nav_res.is_alive(wid)) {
                step_cycles++;
                // stepper_state_t &step_state = stepper_states.at(i);

                if constexpr (std::same_as<actor_chain_t, actor_chain<>>) {
                    // Advance track and transport Jacobian
                    take_step<stepper_t, actor_chain_t>(
                        cfg.stepping, nav_res, step_state,
                        actor_state_refs[stage]);

                    DETRAY_INFO_DEVICE(
                        "STEP %ld: Thread %d: path length %f mm (stage %d, "
                        "data %d, res %d)",
                        step - stage, block.thread_rank(),
                        step_state.path_length(), stage, i, j);
                } else {
                    // Alternate between stepping and acting
                    if (prop_phase[stage] % 2 == 0) {
                        take_step<stepper_t, actor_chain_t>(
                            cfg.stepping, nav_res, step_state,
                            actor_state_refs[stage]);

                        DETRAY_INFO_DEVICE(
                            "STEP %ld: Thread %d: path length %f mm (stage "
                            "%d, "
                            "data %d, res %d)",
                            step - stage, block.thread_rank(),
                            step_state.path_length(), stage, i, j);
                    } else {
                        run_actors<stepper_t, actor_chain_t>(
                            nav_res, step_state, actor_state_refs[stage]);

                        // Do some dummy work
                        const auto trk = step_state();
                        take_step<stepper_t, actor_chain_t>(
                            cfg.stepping, nav_res, step_state,
                            actor_state_refs[stage]);
                        step_state() = trk;

                        DETRAY_INFO_DEVICE(
                            "STEP %ld: Thread %d: Run actors (stage %d, "
                            "data "
                            "%d, res %d)",
                            step - stage, block.thread_rank(), stage, i, j);
                    }
                }

                [[maybe_unused]] const auto &pos = step_state.m_track.pos();
                DETRAY_INFO_DEVICE(
                    "STEP %ld: Thread %d (Setp): pos [%f, %f, %f] [r: %f, "
                    "z: "
                    "%f] "
                    "(stage %d, "
                    "data %d, res %d)",
                    step - stage, block.thread_rank(),
                    *(step_state.m_track.m_vector[0].data()),
                    *(step_state.m_track.m_vector[0].data() + 1),
                    *(step_state.m_track.m_vector[0].data() + 2),
                    vector::perp(pos), pos[2], stage, i, j);
            }

            //__syncwarp();

            //  Publish stepper/actor results (local -> shared memory)
            //  TODO: Make asynchronous by having SoA global state
            // stepper_results[j].free_param(step_state(), wid);
            cuda::memcpy_async(stepper_results[j].pos0.data() + wid,
                               step_state.m_track.m_vector[0].data(),
                               sizeof(scalar_t), propagation);
            cuda::memcpy_async(stepper_results[j].pos1.data() + wid,
                               step_state.m_track.m_vector[0].data() + 1,
                               sizeof(scalar_t), propagation);
            cuda::memcpy_async(stepper_results[j].pos2.data() + wid,
                               step_state.m_track.m_vector[0].data() + 2,
                               sizeof(scalar_t), propagation);
            cuda::memcpy_async(stepper_results[j].dir0.data() + wid,
                               step_state.m_track.m_vector[0].data() + 4,
                               sizeof(scalar_t), propagation);
            cuda::memcpy_async(stepper_results[j].dir1.data() + wid,
                               step_state.m_track.m_vector[0].data() + 5,
                               sizeof(scalar_t), propagation);
            cuda::memcpy_async(stepper_results[j].dir2.data() + wid,
                               step_state.m_track.m_vector[0].data() + 6,
                               sizeof(scalar_t), propagation);

            cuda::std::chrono::high_resolution_clock::time_point cycle_t2 =
                cuda::std::chrono::high_resolution_clock::now();

            // Trigger navigation
            propagation.producer_commit();

            // Flip the propagation phase for this track (stepping <->
            // actors)
            prop_phase[stage] = (prop_phase[stage] + 1) % 2;
            // Flip stage and run on the other track
            stage = (stage + 1) % n_stages;

            const auto cycle_time = cuda::std::chrono::duration_cast<
                cuda::std::chrono::duration<double>>(cycle_t2 - cycle_t1);
            step_work += cycle_time.count();
        }

        DETRAY_INFO_DEVICE("Thread %d, block %d (Step): exit pipeline",
                           block.thread_rank(), blockIdx.x);

        cuda::std::chrono::high_resolution_clock::time_point step_t2 =
            cuda::std::chrono::high_resolution_clock::now();

        const auto step_time = cuda::std::chrono::duration_cast<
            cuda::std::chrono::duration<double>>(step_t2 - step_t1);
        step_total += step_time.count();
    }

    cuda::std::chrono::high_resolution_clock::time_point prop_t2 =
        cuda::std::chrono::high_resolution_clock::now();

    const auto prop_time =
        cuda::std::chrono::duration_cast<cuda::std::chrono::duration<double>>(
            prop_t2 - prop_t1);

    block.sync();
    cuda::std::chrono::high_resolution_clock::time_point block_t2 =
        cuda::std::chrono::high_resolution_clock::now();
    const auto block_time =
        cuda::std::chrono::duration_cast<cuda::std::chrono::duration<double>>(
            block_t2 - block_t1);

    if (wid == 0) {
        printf("%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", gid,
               blockIdx.x, wpid, nav_cycles, step_cycles, nav_work * 1000.,
               step_work * 1000., nav_mem * 1000., step_mem * 1000.,
               nav_wait * 1000., step_wait * 1000., nav_total * 1000.,
               step_total * 1000., init_time.count() * 1000.,
               prop_time.count() * 1000., block_time.count() * 1000.);
    }
}

template <typename device_detector_t>
device_detector_t *setup_device_detector(
    typename device_detector_t::view_type det_view) {

    // Build a device detector type (the interal pointers and capacities
    // refer to the already allocated vecmem device buffers)
    device_detector_t device_det{det_view};

    // Allocate global memory space for the device detector to be shared
    // by kernels
    device_detector_t *pinned_detector_ptr{nullptr};

    DETRAY_CUDA_ERROR_CHECK(cudaHostAlloc((void **)&pinned_detector_ptr,
                                          sizeof(device_detector_t),
                                          cudaHostAllocPortable));

    DETRAY_CUDA_ERROR_CHECK(cudaMemcpy(pinned_detector_ptr, &device_det,
                                       sizeof(device_detector_t),
                                       cudaMemcpyHostToDevice));

    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    return pinned_detector_ptr;
}

template <typename device_detector_t>
void release_device_detector(device_detector_t *pinned_detector_ptr) {
    DETRAY_CUDA_ERROR_CHECK(cudaFreeHost(pinned_detector_ptr));
}

template <typename actor_chain_t>
typename actor_chain_t::state_tuple *setup_actor_states(
    typename actor_chain_t::state_tuple *input_actor_states) {

    // Copy the actor state blueprint to the device
    using actor_state_t = typename actor_chain_t::state_tuple;
    actor_state_t *pinned_actor_state_ptr{nullptr};

    DETRAY_CUDA_ERROR_CHECK(cudaHostAlloc((void **)&pinned_actor_state_ptr,
                                          sizeof(actor_state_t),
                                          cudaHostAllocPortable));
    DETRAY_CUDA_ERROR_CHECK(
        cudaMemcpy(pinned_actor_state_ptr, input_actor_states,
                   sizeof(actor_state_t), cudaMemcpyHostToDevice));

    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    return pinned_actor_state_ptr;
}

template <typename actor_chain_t>
void release_actor_states(
    typename actor_chain_t::state_tuple *pinned_actor_state_ptr) {
    DETRAY_CUDA_ERROR_CHECK(cudaFreeHost(pinned_actor_state_ptr));
}

template <typename navigator_t, typename stepper_t, typename actor_chain_t>
void run_propagation_kernel(
    const propagation::config &cfg,
    const typename navigator_t::detector_type *pinned_detector_ptr,
    typename stepper_t::magnetic_field_type field_view,
    vecmem::data::vector_view<typename stepper_t::free_track_parameters_type>
        tracks_view,
    vecmem::data::vector_view<typename navigator_t::state>
        navigator_states_view,
    vecmem::data::vector_view<typename stepper_t::state> stepper_states_view,
    typename actor_chain_t::state_tuple *pinned_actor_state_ptr) {

    int thread_dim = 256;
    int block_dim =
        (static_cast<int>(tracks_view.size()) + thread_dim - 1) / thread_dim;

    DETRAY_INFO_HOST_DEVICE("# Tracks: %ld", tracks_view.size());
    DETRAY_INFO_HOST_DEVICE("# threads per block: %d", thread_dim);
    DETRAY_INFO_HOST_DEVICE("# blocks: %d", block_dim);
    DETRAY_INFO_HOST_DEVICE("# threads: %d", thread_dim * block_dim);

    std::chrono::high_resolution_clock::time_point t1 =
        std::chrono::high_resolution_clock::now();

    // run the propagation loop
    propagation_kernel<navigator_t, stepper_t, actor_chain_t>
        <<<block_dim, thread_dim>>>(
            cfg, pinned_detector_ptr, field_view, tracks_view,
            navigator_states_view, stepper_states_view, pinned_actor_state_ptr);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    std::chrono::high_resolution_clock::time_point t2 =
        std::chrono::high_resolution_clock::now();

    const auto total_time =
        std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    std::cout << "Propagation kernel took: " << total_time.count() * 1000.
              << "ms\n";
}

/// Macro declaring the template instantiations for the different detector types
#define DECLARE_PROPAGATOR(METADATA, CHAIN, FIELD)                       \
                                                                         \
    template void run_propagation_kernel<                                \
        navigator_type<METADATA>, stepper_type<METADATA, FIELD>,         \
        CHAIN<detector<METADATA>::algebra_type>>(                        \
        const propagation::config &,                                     \
        const detector<METADATA, device_container_types> *,              \
        covfie::field_view<FIELD>,                                       \
        vecmem::data::vector_view<                                       \
            free_track_parameters<detector<METADATA>::algebra_type>>,    \
        vecmem::data::vector_view<navigator_type<METADATA>::state>,      \
        vecmem::data::vector_view<stepper_type<METADATA, FIELD>::state>, \
        typename CHAIN<detector<METADATA>::algebra_type>::state_tuple *);

#define DECLARE_ACTOR_CHAIN_SETUP(METADATA, CHAIN)                           \
                                                                             \
    template typename CHAIN<detector<METADATA>::algebra_type>::state_tuple * \
    setup_actor_states<CHAIN<detector<METADATA>::algebra_type>>(             \
        typename CHAIN<detector<METADATA>::algebra_type>::state_tuple *);    \
                                                                             \
    template void                                                            \
    release_actor_states<CHAIN<detector<METADATA>::algebra_type>>(           \
        typename CHAIN<detector<METADATA>::algebra_type>::state_tuple *);

#define DECLARE_DETECTOR_ALLOCATION(METADATA)                                \
                                                                             \
    template detector<METADATA, device_container_types>                      \
        *setup_device_detector<detector<METADATA, device_container_types>>(  \
            typename detector<METADATA, device_container_types>::view_type); \
                                                                             \
    template void                                                            \
    release_device_detector<detector<METADATA, device_container_types>>(     \
        detector<METADATA, device_container_types> *);

DECLARE_PROPAGATOR(benchmarks::default_metadata, empty_chain, const_field_t)
DECLARE_PROPAGATOR(benchmarks::default_metadata, default_chain, const_field_t)

DECLARE_PROPAGATOR(benchmarks::toy_metadata, empty_chain, const_field_t)
DECLARE_PROPAGATOR(benchmarks::toy_metadata, default_chain, const_field_t)

// Declare only once per algebra type
DECLARE_ACTOR_CHAIN_SETUP(benchmarks::toy_metadata, empty_chain)
DECLARE_ACTOR_CHAIN_SETUP(benchmarks::toy_metadata, default_chain)

DECLARE_DETECTOR_ALLOCATION(benchmarks::default_metadata)
DECLARE_DETECTOR_ALLOCATION(benchmarks::toy_metadata)

}  // namespace detray
