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
    const unsigned int wp_rank{block.thread_rank() % 32u};

    assert(wp_rank < 32u);

    constexpr navigator_t navigator{};

    navigator.update(stepper_result.free_param(wp_rank), navigation, cfg, gtx);
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
    const unsigned int wp_rank{block.thread_rank() % 32u};

    assert(wp_rank < 32u);

    DETRAY_INFO_DEVICE("Warp rank in 'take_step': %d\nNav result dist: %f",
                       wp_rank, nav_result(wp_rank));

    constexpr stepper_t stepper{};

    // Set access to the volume material for the stepper
    /*auto vol = nav_result.current_volume(wp_rank);
    const material<scalar_t> *vol_mat_ptr =
        vol.has_material() ? vol.material_parameters(step_state().pos())
                           : nullptr;*/

    // Break automatic step size scaling by the stepper when a
    // surface was reached and whenever the navigation is
    // (re-)initialized
    const bool reset_stepsize{nav_result.is_on_surface(wp_rank) || true};

    // Take the step
    stepper.step(nav_result(wp_rank), step_state, cfg, reset_stepsize, nullptr);

    // Reduce navigation trust level according to stepper update
    // typename stepper_t::policy_type{}(stepping.policy_state(),
    //                                    propagation);
    nav_result.set_high_trust(wp_rank);
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
__global__ void __launch_bounds__(256, 2) propagation_kernel(
    const propagation::config cfg,
    const typename navigator_t::detector_type *pinned_detector_ptr,
    typename stepper_t::magnetic_field_type field_view,
    vecmem::data::vector_view<typename stepper_t::free_track_parameters_type>
        tracks_view,
    vecmem::data::vector_view<typename navigator_t::state>
        navigator_states_view,
    vecmem::data::vector_view<stepping::result<
        typename navigator_t::detector_type::algebra_type, 32u>>
        stepper_res_view,
    typename actor_chain_t::state_tuple *pinned_actor_state_ptr) {

    using detector_device_t = typename navigator_t::detector_type;
    using metadata_t = typename detector_device_t::metadata;
    using algebra_t = typename detector_device_t::algebra_type;

    using track_t = typename stepper_t::free_track_parameters_type;
    using stepper_state_t = typename stepper_t::state;
    using navigation_state_t = typename navigator_t::state;

    assert(blockDim.y == blockDim.z == 1);

    constexpr int n_propagation_steps{25};
    constexpr unsigned int n_warps{256 / 32u};
    auto block = cooperative_groups::this_thread_block();

    const unsigned int gid = threadIdx.x + blockIdx.x * blockDim.x;
    const unsigned int bid = block.thread_rank();

    // Divide the block into warp pairs, where work is divided between warps
    // into nav <-> step/act
    assert(block.size() % 64 == 0);
    cooperative_groups::thread_block_tile<64> warp_pair =
        cooperative_groups::tiled_partition<64>(block);

    /// One navigation and one stepping/acting thread work on 2 tracks (2:2)
    const vecmem::device_vector<track_t> tracks(tracks_view);

    /// The navigation and stepper result types
    __shared__ navigation::result<metadata_t, 32u> navigator_results[n_warps];
    __shared__ stepping::result<algebra_t, 32u> stepper_results[n_warps];

    // Create a pipeline.
    constexpr auto scope{cuda::thread_scope_block};
    constexpr int n_stages{2};

    // One pipeline per warp pair: 256 / 64 = 4
    __shared__ cuda::pipeline_shared_state<scope, n_stages>
        shared_states[n_warps / 2u];

    // Get the index of the warp pair in the block to assign the shared state
    const unsigned int wid{bid / 64u};
    DETRAY_INFO_DEVICE("Thread %d (%d) in block %d: warp pair index: %d", gid,
                       bid, blockIdx.x, wid);

    // Create the partitioned propagation pipeline for the warp pair
    const bool does_navigation{warp_pair.thread_rank() < warp_pair.size() / 2u};
    auto thread_role = does_navigation ? cuda::pipeline_role::consumer
                                       : cuda::pipeline_role::producer;
    auto propagation =
        cuda::make_pipeline(warp_pair, &shared_states[wid], thread_role);
    constexpr int n_pipline_cycles{4 * n_propagation_steps};

    // Map the id in the warp back onto the same data offset
    const unsigned int data_idx{does_navigation ? gid : gid - 32u};

    // Have both warps compute one of the two initial navigation runs
    const typename detector_device_t::geometry_context gctx{};
    // Setup the navigator and stepper results buffer
    vecmem::device_vector<navigation_state_t> navigator_states{
        navigator_states_view};
    assert(navigator_states.size() == tracks.size());
    navigation_state_t tmp_nav_state{*pinned_detector_ptr};

    {
        // Navigation stage: which of the two tracks to navigate
        const unsigned int stage{does_navigation ? 0u : 1u};
        const unsigned int wp_rank{warp_pair.thread_rank() - stage * 32u};
        const unsigned int i{data_idx + stage * 32u};
        const unsigned int j{i / 32u - blockIdx.x * n_warps};

        if (i < tracks.size()) {
            constexpr navigator_t navigator{};

            navigation_state_t &nav_state_ref =
                does_navigation ? tmp_nav_state : navigator_states[gid];
            nav_state_ref = navigation_state_t{*pinned_detector_ptr};
            navigator.init(tracks.at(i), nav_state_ref, cfg.navigation, gctx);

            if (wp_rank == 0) {
                navigator_results[j] = {};
            }
            navigator_results[j].m_volume_index[wp_rank] =
                nav_state_ref.volume();
            navigator_results[j].m_dist_to_next[wp_rank] = nav_state_ref();
            navigator_results[j].m_status[wp_rank] = nav_state_ref.status();

            DETRAY_INFO_DEVICE(
                "Thread %d: Navigation init: %f (stage %d, data %d)",
                block.thread_rank(), nav_state_ref(), stage, i);

            assert(nav_state_ref.trust_level() ==
                   navigation::trust_level::e_full);
            assert(nav_state_ref.is_alive());
        }
    }
    warp_pair.sync();

    ///
    /// Run propagation steps
    ///

    // The first half of the warp pair performs navigation, the second half
    // stepping and actor calls
    if (does_navigation) {

        if (data_idx + 32u < tracks.size()) {
            DETRAY_INFO_DEVICE(
                "Thread %d (%d) in block %d is doing navigation (tracks: %d, "
                "%d)",
                gid, bid, blockIdx.x, data_idx, data_idx + 32u);
        } else if (data_idx < tracks.size()) {
            DETRAY_INFO_DEVICE(
                "Thread %d (%d) in block %d is doing navigation (tracks: %d, "
                "-)",
                gid, bid, blockIdx.x, data_idx);
        }

        // Fetch the navigation states into local memory for coalesced access
        std::array<navigation_state_t, 2u> navigation{
            tmp_nav_state, navigator_states[gid + 32u]};

        // Propagation loop
        int stage = 0;  //< which of the two tracks to navigate
        const unsigned int wp_rank{warp_pair.thread_rank()};
        for (std::size_t step = 0; step < n_pipline_cycles; ++step) {
            const unsigned int i{data_idx + stage * 32u};
            const unsigned int j{i / 32u - blockIdx.x * n_warps};

            // Wait for the stepper/actor to commit next step/action to pipeline
            // __syncwarp();
            propagation.consumer_wait();

            auto &nav_res = navigator_results[j];

            if (i < tracks.size() && nav_res.is_alive(wp_rank)) {
                if (nav_res.trust_level(wp_rank) !=
                    navigation::trust_level::e_full) {
                    navigation[stage].set_high_trust();
                }
                navigate<navigator_t>(cfg.navigation, gctx, stepper_results[j],
                                      navigation[stage]);

                DETRAY_INFO_DEVICE(
                    "STEP %ld: Thread %d: Dist %f (stage %d, data %d)",
                    step - stage, block.thread_rank(), navigation[stage](),
                    stage, i);

                nav_res.m_volume_index[wp_rank] = navigation[stage].volume();
                nav_res.m_dist_to_next[wp_rank] = navigation[stage]();
                nav_res.m_status[wp_rank] = navigation[stage].status();
            }

            // Publish navigator results
            // __syncwarp();
            propagation.consumer_release();

            // Flip stage and navigate the other track
            stage = (stage + 1) % n_stages;
        }

        DETRAY_INFO_HOST_DEVICE(
            "Thread %d, block %d (Nav): Navigation finished",
            block.thread_rank(), blockIdx.x);
    } else {

        if (data_idx + 32u < tracks.size()) {
            DETRAY_INFO_DEVICE(
                "Thread %d (%d) in block %d is doing stepping (tracks: "
                "%d, "
                "%d)",
                gid, bid, blockIdx.x, data_idx, data_idx + 32u);
        } else if (data_idx < tracks.size()) {
            DETRAY_INFO_DEVICE(
                "Thread %d (%d) in block %d is doing stepping (tracks: "
                "%d, "
                "-)",
                gid, bid, blockIdx.x, data_idx);
        }

        const unsigned int wp_rank{warp_pair.thread_rank() - 32u};

        // Initial track position and direction (actors might change it!)
        std::array<stepper_state_t, 2u> stepping{stepper_state_t{field_view},
                                                 stepper_state_t{field_view}};

        // Make a copy of the initial actor states (blueprint) into local memory
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
            propagation.producer_acquire();

            if (i < tracks.size()) {
                auto &nav_res = navigator_results[j];
                stepping[stage]() = tracks.at(data_idx);

                // Check actor update
                run_actors<stepper_t, actor_chain_t>(nav_res, stepping[stage],
                                                     actor_state_refs[stage]);

                take_step<stepper_t, actor_chain_t>(cfg.stepping, nav_res,
                                                    stepping[stage],
                                                    actor_state_refs[stage]);

                DETRAY_INFO_DEVICE(
                    "Thread %d: Fill pipeline: Run actors (stage %d, data "
                    "%d)",
                    block.thread_rank(), stage, i);

                // Publish stepper results
                stepper_results[j].free_param(stepping[stage](), wp_rank);
            }

            // __syncwarp();
            propagation.producer_commit();
        }

        // Propagation loop
        int stage = 0;                      //< which track to advance
        std::array<int, 2u> prop_phase{0};  //< whether to run actors or stepper
        for (std::size_t step = 0; step < n_pipline_cycles; ++step) {
            const unsigned int i{data_idx + stage * 32u};
            const unsigned int j{i / 32u - blockIdx.x * n_warps};

            // Wait for navigation to finish
            // __syncwarp();
            propagation.producer_acquire();
            auto &nav_res = navigator_results[j];

            if (i < tracks.size() && nav_res.is_alive(wp_rank)) {
                if (prop_phase[stage] % 2 == 0) {
                    take_step<stepper_t, actor_chain_t>(
                        cfg.stepping, nav_res, stepping[stage],
                        actor_state_refs[stage]);

                    DETRAY_INFO_DEVICE(
                        "STEP %ld: Thread %d: path length %f (stage %d, data "
                        "%d)",
                        step - stage, block.thread_rank(),
                        stepping[stage].path_length(), stage, i);

                    // Publish stepper results
                    stepper_results[j].free_param(stepping[stage](), wp_rank);
                } else {
                    // Check actor update
                    run_actors<stepper_t, actor_chain_t>(
                        nav_res, stepping[stage], actor_state_refs[stage]);

                    take_step<stepper_t, actor_chain_t>(
                        cfg.stepping, nav_res, stepping[stage],
                        actor_state_refs[stage]);

                    DETRAY_INFO_DEVICE(
                        "STEP %ld: Thread %d: Run actors (stage %d, data "
                        "%d)",
                        step - stage, block.thread_rank(), stage, i);
                }
            }

            // Trigger navigation
            // __syncwarp();
            propagation.producer_commit();

            // Flip the propagation phase for this track (stepping <-> actors)
            prop_phase[stage] = (prop_phase[stage] + 1) % 2;
            // Flip stage and run on the other track
            stage = (stage + 1) % n_stages;
        }

        DETRAY_INFO_HOST_DEVICE("Thread %d, block %d (Step): Stepping finished",
                                block.thread_rank(), blockIdx.x);
    }
}

template <typename device_detector_t>
device_detector_t *setup_device_detector(
    typename device_detector_t::view_type det_view) {

    // Build a device detector type (the interal pointers and capacities
    // refer to the already allocated vecmem device buffers)
    device_detector_t device_det{det_view};

    // Allocate global memory space for the device detector to be shared by
    // kernels
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
    vecmem::data::vector_view<stepping::result<
        typename navigator_t::detector_type::algebra_type, 32u>>
        stepper_res_view,
    typename actor_chain_t::state_tuple *pinned_actor_state_ptr) {

    int thread_dim = 256;
    int block_dim =
        (static_cast<int>(tracks_view.size()) + thread_dim - 1) / thread_dim;

    DETRAY_INFO_HOST_DEVICE("# Tracks: %ld", tracks_view.size());
    DETRAY_INFO_HOST_DEVICE("# threads per block: %d", thread_dim);
    DETRAY_INFO_HOST_DEVICE("# blocks: %d", block_dim);
    DETRAY_INFO_HOST_DEVICE("# threads: %d", thread_dim * block_dim);

    // run the propagation loop
    propagation_kernel<navigator_t, stepper_t, actor_chain_t>
        <<<block_dim, thread_dim>>>(cfg, pinned_detector_ptr, field_view,
                                    tracks_view, navigator_states_view,
                                    stepper_res_view, pinned_actor_state_ptr);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/// Macro declaring the template instantiations for the different detector types
#define DECLARE_PROPAGATOR(METADATA, CHAIN, FIELD)                            \
                                                                              \
    template void run_propagation_kernel<                                     \
        navigator_type<METADATA>, stepper_type<METADATA, FIELD>,              \
        CHAIN<detector<METADATA>::algebra_type>>(                             \
        const propagation::config &,                                          \
        const detector<METADATA, device_container_types> *,                   \
        covfie::field_view<FIELD>,                                            \
        vecmem::data::vector_view<                                            \
            free_track_parameters<detector<METADATA>::algebra_type>>,         \
        vecmem::data::vector_view<navigator_type<METADATA>::state>,           \
        vecmem::data::vector_view<                                            \
            detray::stepping::result<detector<METADATA>::algebra_type, 32u>>, \
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
