/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/benchmarks/device/cuda/propagator.hpp"
#include "detray/definitions/detail/cuda_definitions.hpp"

// CUDA include(s)
#include <cooperative_groups.h>

#include <cuda/barrier>

namespace detray {

/// Initialize the stepper, navigator and actor states on device
template <typename navigator_t, typename stepper_t, typename actor_chain_t>
__global__ void __launch_bounds__(256, 4) propagation_init_kernel(
    const propagation::config *pinned_cfg,
    typename navigator_t::detector_type::view_type det_view,
    typename stepper_t::magnetic_field_type field_view,
    vecmem::data::vector_view<typename stepper_t::free_track_parameters_type>
        tracks_view,
    vecmem::data::vector_view<typename stepper_t::state> stepper_states_view,
    vecmem::data::vector_view<typename navigator_t::state>
        navigator_states_view,
    vecmem::data::vector_view<typename actor_chain_t::state_tuple>
        actor_states_view,
    typename actor_chain_t::state_tuple *pinned_actor_state_ptr) {

    using metadata_t = typename navigator_t::detector_type::metadata;
    using detector_device_t = detector<metadata_t, device_container_types>;
    using algebra_t = typename detector_device_t::algebra_type;

    using track_t = typename stepper_t::free_track_parameters_type;
    using stepper_state_t = typename stepper_t::state;
    using navigation_state_t = typename navigator_t::state;
    using propagation_state_t =
        typename cuda_propagation<navigator_t, stepper_t, actor_chain_t>::state;

    const detector_device_t det(det_view);
    const vecmem::device_vector<track_t> tracks(tracks_view);
    vecmem::device_vector<stepper_state_t> stepper_states(stepper_states_view);
    vecmem::device_vector<navigation_state_t> navigation_states(
        navigator_states_view);
    vecmem::device_vector<typename actor_chain_t::state_tuple> actor_states(
        actor_states_view);

    constexpr navigator_t navigator{};
    constexpr actor_chain_t run_actors{};

    typename detector_device_t::geometry_context gctx{};

    // Setup some small amount of data in shared memory
    __shared__ navigation::config nav_cfg[1];
    __shared__ typename actor_chain_t::state_tuple actor_states_sh[1];

    auto block = cooperative_groups::this_thread_block();
    if (block.thread_rank() == 0) {
        nav_cfg[0] = pinned_cfg->navigation;
    }

    if (block.thread_rank() == 0) {
        actor_states_sh[0] = *pinned_actor_state_ptr;
    }
    block.sync();

    /// Every thread initializes one set of states
    int gid = threadIdx.x + blockIdx.x * blockDim.x;
    if (gid >= tracks.size()) {
        return;
    }

    // Create the stepper state:
    // The track gets copied into the stepper state, so that the
    // original track sample vector remains unchanged
    stepper_state_t &stepping = stepper_states.at(gid);
    stepping = stepper_state_t{tracks.at(gid), field_view};

    // Create navigator state
    navigation_state_t &navigation = navigation_states.at(gid);
    navigation = navigation_state_t{det};

    // Temporary propagation state
    propagation_state_t propagation{stepping, navigation, gctx};

    // Create the actor states on a fresh copy
    auto actor_states_ref =
        actor_chain_t::setup_actor_states(actor_states.at(gid));
    actor_states_ref = actor_states_sh[0];

    // Initialize the navigation
    navigator.init(stepping(), navigation, *nav_cfg, gctx);
    assert(navigation.is_alive());

    // Run all registered actors/aborters
    run_actors(actor_states_ref, propagation);
    assert(!stepping().is_invalid());
    assert(!stepping.bound_params().is_invalid());

    // Update the navigation information, in case the actors changed the track
    navigator.update(stepping(), navigation, *nav_cfg, gctx);
    assert(navigation.is_alive());
}

template <typename collection_t>
__device__ std::size_t scan(collection_t &c) {

    for (std::size_t i = 0; i < c.size(); ++i) {
    }
}

/// Run the stepper
template <typename navigator_t, typename stepper_t, typename actor_chain_t>
__device__ inline void take_step(
    const stepping::config *cfg, cooperative_groups::thread_block block,
    const std::array<int, 2> queues[],
    ::cuda::barrier<::cuda::thread_scope_block> needs_stepping[],
    ::cuda::barrier<::cuda::thread_scope_block> needs_navigation[],
    vecmem::device_vector<typename stepper_t::state> stepper_states,
    // vecmem::device_vector<unsigned int> stepper_res_view,
    vecmem::device_vector<typename navigator_t::state> navigation_states
    // vecmem::device_vector<unsigned int> navigator_res_view
) {

    using barrier_t = ::cuda::barrier<::cuda::thread_scope_block>;

    using detector_t = typename navigator_t::detector_type;
    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;

    using stepper_state_t = typename stepper_t::state;
    using navigation_state_t = typename navigator_t::state;
    using propagation_state_t =
        typename cuda_propagation<navigator_t, stepper_t, actor_chain_t>::state;

    constexpr stepper_t stepper{};
    constexpr bool reset_stepsize{true};
    typename detector_t::geometry_context gctx{};

    assert(cfg != nullptr);
    assert(stepper_states.size() == navigation_states.size());
    assert(static_cast<int>(stepper_states.size()) > 512 * blockIdx.x);

    // Block-local thread index
    const int thread_idx = block.thread_rank();

    // Does this thread have work in either queue ?
    const int state_idx_0{thread_idx + queues[0][0]};
    const int state_idx_1{thread_idx + queues[1][0]};

    const bool has_queue_0{state_idx_1 < queues[0][1]};
    const bool has_queue_1{state_idx_1 < queues[1][1]};

    // Run while either at least one of the two navigation streams monitored by
    // this thread is alive
    int i{0};
    while ((has_queue_0 && navigation_states.at(state_idx_0).is_alive() ||
            (has_queue_1 && navigation_states.at(state_idx_1).is_alive())) &&
           i < 1000) {
        // Queue for this iteration
        const int queue_idx{i % 2};
        const std::array<int, 2> &current_queue{queues[queue_idx]};
        // Which state to get
        auto state_idx = thread_idx + current_queue[0];

        printf("Thread %d (Step): Wait for navigation (queue: %d)\n",
               thread_idx, queue_idx);
        needs_stepping[queue_idx].arrive_and_wait();

        printf("Thread %d (Step): Run state index %d (queue: %d)\n", thread_idx,
               state_idx, queue_idx);
        if (state_idx < current_queue[1]) {
            navigation_state_t &navigation = navigation_states.at(state_idx);

            if (navigation.is_alive()) {
                printf("Thread %d (Step): Do stepping (queue: %d)\n",
                       thread_idx, queue_idx);

                stepper_state_t &stepping = stepper_states.at(state_idx);
                const auto &track = stepping();

                /*const auto vol = navigation.get_volume();
                const material<scalar_t> *vol_mat_ptr =
                    vol.has_material() ? vol.material_parameters(track.pos())
                                       : nullptr;*/

                stepper.step(navigation(), stepping, *cfg, reset_stepsize,
                             /*vol_mat_ptr*/ nullptr);

                propagation_state_t propagation{stepping, navigation, gctx};
                typename stepper_t::policy_type{}(stepping.policy_state(),
                                                  propagation);
            }
        }

        if (!(has_queue_0 && navigation_states.at(state_idx_0).is_alive()) &&
            !(has_queue_1 && navigation_states.at(state_idx_1).is_alive())) {
            printf("Thread %d (Step): Drop out of propagation\n", thread_idx);
            needs_navigation[0].arrive_and_drop();
            needs_navigation[1].arrive_and_drop();

            needs_stepping[0].arrive_and_wait();
            needs_stepping[1].arrive_and_wait();
        } else {
            printf("Thread %d (Step): Needs navigation (queue: %d)\n",
                   thread_idx, queue_idx);
            barrier_t::arrival_token token =
                needs_navigation[queue_idx].arrive();
        }
        ++i;
        if (blockIdx.x == 1) {
            printf("Thread %d, block %d (Step): Iteration %d\n", thread_idx,
                   blockIdx.x, i);
        }
    }
    if (i >= 1000) {
        printf("Thread %d, block %d  (Step): Reached hard limit\n", thread_idx,
               blockIdx.x);
        needs_navigation[0].arrive_and_drop();
        needs_navigation[1].arrive_and_drop();

        needs_stepping[0].arrive_and_wait();
        needs_stepping[1].arrive_and_wait();
    }
    printf("Thread %d, block %d (Step): Stepping finished\n", thread_idx,
           blockIdx.x);
}

/// Run the navigator
template <typename navigator_t, typename stepper_t, typename actor_chain_t>
__device__ void navigate(
    const navigation::config *cfg, cooperative_groups::thread_block block,
    const std::array<int, 2> queues[],
    ::cuda::barrier<::cuda::thread_scope_block> needs_stepping[],
    ::cuda::barrier<::cuda::thread_scope_block> needs_navigation[],
    const vecmem::device_vector<typename stepper_t::state> stepper_states,
    // vecmem::device_vector<unsigned int> stepper_res_view,
    vecmem::device_vector<typename navigator_t::state> navigation_states
    // vecmem::device_vector<unsigned int> navigator_res_view
) {

    using barrier_t = ::cuda::barrier<::cuda::thread_scope_block>;

    using detector_t = typename navigator_t::detector_type;
    using stepper_state_t = typename stepper_t::state;
    using navigation_state_t = typename navigator_t::state;

    assert(cfg != nullptr);
    assert(stepper_states.size() == navigation_states.size());
    assert(static_cast<int>(stepper_states.size()) > 512 * blockIdx.x);

    typename detector_t::geometry_context gctx{};
    constexpr navigator_t navigator{};

    // Block-local thread index (adjusted for the start of the warp
    // specialization)
    const int thread_idx = block.thread_rank() - blockDim.x / 2;
    assert(thread_idx >= 0);

    // Wait for all threads in the block to arrive for stepping
    printf("Thread %d (Nav): Needs stepping (queue: 0)\n", thread_idx);
    barrier_t::arrival_token token1 = needs_stepping[0].arrive();
    printf("Thread %d (Nav): Needs stepping (queue: 1)\n", thread_idx);
    barrier_t::arrival_token token2 = needs_stepping[1].arrive();

    // Does this thread have work in either queue ?
    const int state_idx_0{thread_idx + queues[0][0]};
    const int state_idx_1{thread_idx + queues[1][0]};

    const bool has_queue_0{state_idx_1 < queues[0][1]};
    const bool has_queue_1{state_idx_1 < queues[1][1]};

    if (blockIdx.x == 1) {
        printf(
            "Thread %d, block %d (Nav): has queue 0 %d (state idx %d [%d, "
            "%d])\n",
            thread_idx, blockIdx.x, has_queue_0, state_idx_0, queues[0][0],
            queues[0][1]);
        printf(
            "Thread %d, block %d (Nav): has queue 1 %d (state idx %d [%d, "
            "%d])\n",
            thread_idx, blockIdx.x, has_queue_1, state_idx_1, queues[1][0],
            queues[1][1]);
    }

    // Run while either at least one of the two navigation streams monitored by
    // this thread is alive
    int i{0};
    while ((has_queue_0 && navigation_states.at(state_idx_0).is_alive() ||
            (has_queue_1 && navigation_states.at(state_idx_1).is_alive())) &&
           i < 1000) {
        // Queue for this iteration
        const int queue_idx{i % 2};
        const std::array<int, 2> &current_queue{queues[queue_idx]};
        // Which state to get
        auto state_idx = thread_idx + current_queue[0];

        printf("Thread %d (Nav): Wait for stepping (queue: %d)\n", thread_idx,
               queue_idx);
        // Wait for stepping to fill the queue
        needs_navigation[queue_idx].arrive_and_wait();
        printf("Thread %d (Nav): Run state index %d (queue: %d)\n", thread_idx,
               state_idx, queue_idx);

        // If eligible to do work: navigate
        if (state_idx < current_queue[1]) {

            navigation_state_t &navigation = navigation_states.at(state_idx);

            if (navigation.is_alive()) {
                const auto &stepper_state = stepper_states.at(state_idx);
                navigator.update(stepper_state(), navigation, *cfg, gctx);
            }
        }

        if (!(has_queue_0 && navigation_states.at(state_idx_0).is_alive()) &&
            !(has_queue_1 && navigation_states.at(state_idx_1).is_alive())) {
            printf("Thread %d (Nav): Drop out of propagation\n", thread_idx);
            // No more stepping needed
            needs_stepping[0].arrive_and_drop();
            needs_stepping[1].arrive_and_drop();

            // Wait for stepping to signal 'no more navigation needed'
            needs_navigation[0].arrive_and_wait();
            needs_navigation[1].arrive_and_wait();
        } else {
            printf("Thread %d (Nav): Needs stepping (queue: %d)\n", thread_idx,
                   queue_idx);
            barrier_t::arrival_token token = needs_stepping[queue_idx].arrive();
        }
        ++i;
        if (blockIdx.x == 1) {
            printf("Thread %d, block %d (Nav): Iteration %d\n", thread_idx,
                   blockIdx.x, i);
        }
    }
    if (i >= 1000) {
        printf("Thread %d, block %d  (Nav): Reached hard limit\n", thread_idx,
               blockIdx.x);
        needs_stepping[0].arrive_and_drop();
        needs_stepping[1].arrive_and_drop();

        needs_navigation[0].arrive_and_wait();
        needs_navigation[1].arrive_and_wait();
    }

    printf("Thread %d, block %d (Nav): Navigation finished\n", thread_idx,
           blockIdx.x);
}

/// Specialize the warps to run stepping, navigation and actors independently
template <typename navigator_t, typename stepper_t, typename actor_chain_t>
__global__ void __launch_bounds__(256, 4) propagation_kernel(
    const propagation::config *pinned_cfg,
    vecmem::data::vector_view<typename stepper_t::state> stepper_states_view,
    // vecmem::data::vector_view<unsigned int> stepper_res_view,
    vecmem::data::vector_view<typename navigator_t::state>
        navigator_states_view,
    // vecmem::data::vector_view<unsigned int> navigator_res_view,
    vecmem::data::vector_view<typename actor_chain_t::state_tuple>
        actor_states_view) {

    using barrier_t = ::cuda::barrier<::cuda::thread_scope_block>;
    using algebra_t = typename navigator_t::detector_type::algebra_type;

    /// Register the actor types
    const actor_chain_t run_actors{};

    vecmem::device_vector<typename stepper_t::state> stepper_states(
        stepper_states_view);
    // vecmem::device_vector<unsigned int> stepper_results(stepper_res_view);
    vecmem::device_vector<typename navigator_t::state> navigation_states(
        navigator_states_view);

    int gid = threadIdx.x + blockIdx.x * blockDim.x;

    assert(blockDim.y == blockDim.z == 1);

    __shared__ stepping::config step_cfg[1];
    __shared__ navigation::config nav_cfg[1];

    // In each iteration, handles only half of the tracks (two queues: try to
    // avoid stalls)
    constexpr int queue_length{256};

    // Map the two queues onto the state containers
    using queue_t = std::array<int, 2>;

    __shared__ queue_t queues[2];

    auto block = cooperative_groups::this_thread_block();
    if (block.thread_rank() == 0) {

        printf("This is block %d\n", blockIdx.x);
        step_cfg[0] = pinned_cfg->stepping;
        nav_cfg[0] = pinned_cfg->navigation;

        // Offset for this block's data into the state containers
        const int trk_offset{512 * static_cast<int>(blockIdx.x)};
        // Number of tracks this block handles (comes in chunks of 512 tracks)
        const int n_tracks{math::min(
            static_cast<int>(stepper_states.size()) - trk_offset, 512)};

        const queue_t first_queue = {
            trk_offset, trk_offset + math::min(n_tracks, queue_length)};
        const queue_t second_queue = {
            trk_offset + queue_length,
            trk_offset + queue_length +
                math::min(math::max(0, n_tracks - queue_length), queue_length)};

        printf("Block %d queue 0: [%d, %d]\n", blockIdx.x, first_queue[0],
               first_queue[1]);
        printf("Block %d queue 1: [%d, %d]\n", blockIdx.x, second_queue[0],
               second_queue[1]);

        queues[0] = first_queue;
        queues[1] = second_queue;
    }
    block.sync();

    // vecmem::device_vector<unsigned int>
    // navigation_results(navigator_res_view);

    /*if (gid < navigator_states_view.size()) {
        const auto dir = stepper_states.at(gid)().dir();
        printf("%f: [%f, %f, %f]\n", navigation_states.at(gid)(), dir[0],
               dir[1], dir[2]);
    }*/

    // bar[0] and bar[1] track if buffers buffer_0 and buffer_1 are ready to be
    // filled, while bar[2] and bar[3] track if buffers buffer_0 and buffer_1
    // are filled-in respectively
    __shared__ barrier_t bar[4];

    if (block.thread_rank() < 4) {
        init(bar + block.thread_rank(), block.size());
    }
    block.sync();

    const int split_idx = blockDim.x / 2;
    if (block.thread_rank() < split_idx) {
        printf("Thread %d (%d) in block %d doing stepping\n", gid,
               block.thread_rank(), blockIdx.x);
        take_step<navigator_t, stepper_t, actor_chain_t>(
            step_cfg, block, queues, bar, bar + 2, stepper_states,
            navigation_states);
    } else {
        printf("Thread %d (%d) in block %d doing navigation\n", gid,
               block.thread_rank(), blockIdx.x);
        navigate<navigator_t, stepper_t, actor_chain_t>(
            nav_cfg, block, queues, bar, bar + 2, stepper_states,
            navigation_states);
    }
}

typename propagation::config *setup_config(
    const typename propagation::config *input_config) {

    // Copy the config to the device
    propagation::config *pinned_config_ptr{nullptr};

    DETRAY_CUDA_ERROR_CHECK(cudaHostAlloc((void **)&pinned_config_ptr,
                                          sizeof(propagation::config),
                                          cudaHostAllocPortable));

    DETRAY_CUDA_ERROR_CHECK(cudaMemcpy(pinned_config_ptr, input_config,
                                       sizeof(propagation::config),
                                       cudaMemcpyHostToDevice));

    return pinned_config_ptr;
}

void release_config(typename propagation::config *pinned_config_ptr) {
    DETRAY_CUDA_ERROR_CHECK(cudaFreeHost(pinned_config_ptr));
}

template <typename device_detector_t>
typename device_detector_t *allocate_device_detector() {

    // Allocate global memory space for the device detector to be shared by
    // kernels
    device_detector_t *device_detector_ptr{nullptr};

    DETRAY_CUDA_ERROR_CHECK(
        cudaAlloc((void **)&device_detector_ptr, sizeof(device_detector_t)));

    return device_detector_ptr;
}

template <typename device_detector_t>
void release_detector(typename device_detector_t *device_detector_ptr) {
    DETRAY_CUDA_ERROR_CHECK(cudaFree(device_detector_ptr));
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

    return pinned_actor_state_ptr;
}

template <typename actor_chain_t>
void release_actor_states(
    typename actor_chain_t::state_tuple *pinned_actor_state_ptr) {
    DETRAY_CUDA_ERROR_CHECK(cudaFreeHost(pinned_actor_state_ptr));
}

template <typename navigator_t, typename stepper_t, typename actor_chain_t>
void run_propagation_init_kernel(
    const propagation::config *pinned_cfg,
    typename navigator_t::detector_type::view_type det_view,
    typename stepper_t::magnetic_field_type field_view,
    vecmem::data::vector_view<typename stepper_t::free_track_parameters_type>
        tracks_view,
    vecmem::data::vector_view<typename stepper_t::state> stepper_states_view,
    vecmem::data::vector_view<typename navigator_t::state>
        navigator_states_view,
    vecmem::data::vector_view<typename actor_chain_t::state_tuple>
        actor_states_view,
    typename actor_chain_t::state_tuple *pinned_actor_state_ptr) {

    constexpr int thread_dim = 256;
    int block_dim = (tracks_view.size() + thread_dim - 1) / thread_dim;

    // run the test kernel
    propagation_init_kernel<navigator_t, stepper_t, actor_chain_t>
        <<<block_dim, thread_dim>>>(
            pinned_cfg, det_view, field_view, tracks_view, stepper_states_view,
            navigator_states_view, actor_states_view, pinned_actor_state_ptr);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

template <typename navigator_t, typename stepper_t, typename actor_chain_t>
void run_propagation_kernel(
    const propagation::config *pinned_cfg,
    vecmem::data::vector_view<typename stepper_t::state> stepper_states_view,
    vecmem::data::vector_view<unsigned int> stepper_res_view,
    vecmem::data::vector_view<typename navigator_t::state>
        navigator_states_view,
    vecmem::data::vector_view<unsigned int> navigator_res_view,
    vecmem::data::vector_view<typename actor_chain_t::state_tuple>
        actor_states_view) {

    constexpr int thread_dim = 256;
    // One block handles 512 tracks
    int block_dim =
        (stepper_states_view.size() + 2 * thread_dim - 1) / (2 * thread_dim);

    std::cout << "# Tracks: " << stepper_states_view.size() << std::endl;
    std::cout << "# threads per block: " << thread_dim
              << "\n# blocks: " << block_dim
              << "\n# threads: " << thread_dim * block_dim << std::endl;

    // run the propagation loop
    propagation_kernel<navigator_t, stepper_t, actor_chain_t>
        <<<block_dim, thread_dim>>>(pinned_cfg, stepper_states_view,
                                    navigator_states_view, actor_states_view);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/// Macro declaring the template instantiations for the different detector types
#define DECLARE_PROPAGATOR(METADATA, CHAIN, FIELD)                            \
                                                                              \
    template void run_propagation_init_kernel<                                \
        navigator_type<METADATA>, stepper_type<METADATA, FIELD>,              \
        CHAIN<detector<METADATA>::algebra_type>>(                             \
        const propagation::config *, typename detector<METADATA>::view_type,  \
        covfie::field_view<FIELD>,                                            \
        vecmem::data::vector_view<                                            \
            free_track_parameters<detector<METADATA>::algebra_type>>,         \
        vecmem::data::vector_view<                                            \
            typename stepper_type<METADATA, FIELD>::state>,                   \
        vecmem::data::vector_view<typename navigator_type<METADATA>::state>,  \
        vecmem::data::vector_view<                                            \
            typename CHAIN<detector<METADATA>::algebra_type>::state_tuple>,   \
        typename CHAIN<detector<METADATA>::algebra_type>::state_tuple *);     \
                                                                              \
    template void run_propagation_kernel<                                     \
        navigator_type<METADATA>, stepper_type<METADATA, FIELD>,              \
        CHAIN<detector<METADATA>::algebra_type>>(                             \
        const propagation::config *,                                          \
        vecmem::data::vector_view<                                            \
            typename stepper_type<METADATA, FIELD>::state>,                   \
        vecmem::data::vector_view<unsigned int>,                              \
        vecmem::data::vector_view<typename navigator_type<METADATA>::state>,  \
        vecmem::data::vector_view<unsigned int>,                              \
        vecmem::data::vector_view<                                            \
            typename CHAIN<detector<METADATA>::algebra_type>::state_tuple>);  \
                                                                              \
    template detector<METADATA, device_container_types> *                     \
    allocate_device_detector<detector<METADATA, device_container_types>>();   \
                                                                              \
    template void                                                             \
    release_device_detector<detector<METADATA, device_container_types>>(      \
        detector<METADATA, device_container_types> *);                        \
                                                                              \
    #define DECLARE_ACTOR_CHAIN_SETUP(METADATA, CHAIN) template               \
        typename CHAIN<detector<METADATA>::algebra_type>::state_tuple *       \
        setup_actor_states<CHAIN<detector<METADATA>::algebra_type>>(          \
            typename CHAIN<detector<METADATA>::algebra_type>::state_tuple *); \
                                                                              \
    template void                                                             \
    release_actor_states<CHAIN<detector<METADATA>::algebra_type>>(            \
        typename CHAIN<detector<METADATA>::algebra_type>::state_tuple *);

DECLARE_PROPAGATOR(benchmarks::default_metadata, empty_chain, const_field_t)
DECLARE_PROPAGATOR(benchmarks::default_metadata, default_chain, const_field_t)

DECLARE_PROPAGATOR(benchmarks::toy_metadata, empty_chain, const_field_t)
DECLARE_PROPAGATOR(benchmarks::toy_metadata, default_chain, const_field_t)

// Declare only once per algebra type
DECLARE_ACTOR_CHAIN_SETUP(benchmarks::toy_metadata, empty_chain)
DECLARE_ACTOR_CHAIN_SETUP(benchmarks::toy_metadata, default_chain)

}  // namespace detray
