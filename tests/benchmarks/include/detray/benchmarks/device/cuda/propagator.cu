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
__global__ void propagation_init_kernel(
    propagation::config cfg,
    typename navigator_t::detector_type::view_type det_view,
    typename stepper_t::magnetic_field_type field_view,
    vecmem::data::vector_view<typename stepper_t::free_track_parameters_type>
        tracks_view,
    vecmem::data::vector_view<typename stepper_t::state> stepper_states_view,
    vecmem::data::vector_view<typename navigator_t::state>
        navigator_states_view,
    vecmem::data::vector_view<typename actor_chain_t::state_tuple>
        actor_states_view,
    typename actor_chain_t::state_tuple *device_actor_state_ptr) {

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
    const unsigned int n_tracks{tracks.size()};

    /// Every thread initializes one set of states
    int gid = threadIdx.x + blockIdx.x * blockDim.x;
    if (gid >= n_tracks) {
        return;
    }

    const auto track = tracks.at(gid);

    // Create the stepper state:
    // The track gets copied into the stepper state, so that the
    // original track sample vector remains unchanged
    stepper_states.at(gid) = stepper_state_t{track, field_view};
    stepper_state_t &stepping = stepper_states.at(gid);

    // Create navigator state
    navigation_states.at(gid) = navigation_state_t{det};
    navigation_state_t &navigation = navigation_states.at(gid);

    // Temporary propagation state
    propagation_state_t propagation{stepping, navigation, gctx};

    // Create the actor states on a fresh copy
    actor_states.at(gid) = *device_actor_state_ptr;
    auto actor_states_ref =
        actor_chain_t::setup_actor_states(actor_states.at(gid));

    // Initialize the navigation
    navigator.init(track, navigation, cfg.navigation, gctx);
    assert(navigation.is_alive());

    // Run all registered actors/aborters after navigation init
    run_actors(actor_states_ref, propagation);
    assert(!stepping().is_invalid());

    // Find next candidate if the actors changed the track parameters
    navigator.update(track, navigation, cfg.navigation, gctx);
    assert(navigation.is_alive());
}

template <typename collection_t>
__device__ std::size_t scan(collection_t &c) {

    for (std::size_t i = 0; i < c.size(); ++i) {
    }
}

/// Run the stepper
template <typename stepper_state_t, typename navigator_state_t>
__device__ inline void do_step(
    ::cuda::barrier<::cuda::thread_scope_block> ready[],
    ::cuda::barrier<::cuda::thread_scope_block> filled[],
    vecmem::device_vector<stepper_state_t> stepper_states_view,
    // vecmem::device_vector<unsigned int> stepper_res_view,
    vecmem::device_vector<navigator_state_t> navigator_states_view,
    // vecmem::device_vector<unsigned int> navigator_res_view,
    int thread_idx) {

    using barrier_t = ::cuda::barrier<::cuda::thread_scope_block>;

    std::size_t n_tracks{stepper_states_view.size()};

    // Wait until navigation is done
    // needs_navigation[0].arrive_and_wait();

    // Set access to the volume material for the stepper
    /*auto vol = navigation.get_volume();
    const material<scalar_type> *vol_mat_ptr =
        vol.has_material() ? vol.material_parameters(track.pos()) : nullptr;

    // Break automatic step size scaling by the stepper when a surface
    // was reached and whenever the navigation is (re-)initialized
    const bool reset_stepsize{navigation.is_on_surface() || is_init};
    // Take the step
    propagation._heartbeat &=
        m_stepper.step(navigation(), stepping, m_cfg.stepping,
                       reset_stepsize, vol_mat_ptr);

    // Reduce navigation trust level according to stepper update
    typename stepper_t::policy_type{}(stepping.policy_state(),
    propagation);*/

    // Ready for navigation (non-blocking: perform stepping for next track)
    // auto token = needs_stepping[0].arrive();

    for (int i = 0; i < (n_tracks / 2); ++i) {

        printf("Thread %d (Step): Wait at barrier %d\n", thread_idx, i % 2);
        ready[i % 2].arrive_and_wait(); /* wait for buffer_(i%2) to be ready to
                                           be filled */
        /* produce, i.e., fill in, buffer_(i%2)  */
        printf("Thread %d (Step): Stepping barrier %d\n", thread_idx, i % 2);
        barrier_t::arrival_token token =
            filled[i % 2].arrive(); /* buffer_(i%2) is filled */
        printf("Thread %d (Step): Navigation state ready %d\n", thread_idx,
               i % 2);
    }
    printf("Thread %d (Step): Stepping finished\n", thread_idx);
}

/// Run the navigator
template <typename stepper_state_t, typename navigator_state_t>
__device__ void navigate(
    ::cuda::barrier<::cuda::thread_scope_block> ready[],
    ::cuda::barrier<::cuda::thread_scope_block> filled[],
    vecmem::device_vector<stepper_state_t> stepper_states_view,
    // vecmem::device_vector<unsigned int> stepper_res_view,
    vecmem::device_vector<navigator_state_t> navigator_states_view,
    // vecmem::device_vector<unsigned int> navigator_res_view,
    int thread_idx) {

    using barrier_t = ::cuda::barrier<::cuda::thread_scope_block>;

    printf("Thread %d (Nav): Steper state ready 0\n", thread_idx);
    barrier_t::arrival_token token1 = ready[0].arrive(); /* buffer_0 is ready
    for initial fill */
    printf("Thread %d (Nav): Steper state ready 1\n", thread_idx);
    barrier_t::arrival_token token2 =
        ready[1].arrive(); /* buffer_1 is ready for initial fill */
    for (int i = 0; i < (stepper_states_view.size() / 2); ++i) {

        printf("Thread %d (Nav): Wait at barrier %d\n", thread_idx, i % 2);
        filled[i % 2]
            .arrive_and_wait(); /* wait for buffer_(i%2) to be filled */
        /* consume buffer_(i%2) */
        printf("Thread %d (Nav): Navigating barrier %d\n", thread_idx, i % 2);
        barrier_t::arrival_token token = ready[i % 2].arrive(); /* buffer_(i%2)
          is ready to be re-filled */
        printf("Thread %d (Nav): Steper state ready %d\n", thread_idx, i % 2);
    }
    printf("Thread %d (Nav): Navigation finished\n", thread_idx);
    // Mark all tracks ready for navigation to kick off producer-consumer loop
    /*auto token1 = needs_navigation[0].arrive();

    std::size_t n_tracks{stepper_states_view.size()};

        printf("Thread %d: Wait(Navigation)\n", thread_idx);
        // Wait for stepping to be complete
        needs_stepping[0].arrive_and_wait();

        // Do the navigation
        printf("Thread %d: Navigating\n", thread_idx);

        // Ready for stepping (non-blocking: perform navigation for next track)
        auto token = needs_navigation[0].arrive();*/
}

/// Specialize the warps to run stepping, navigation and actors independently
template <typename navigator_t, typename stepper_t, typename actor_chain_t>
__global__ void propagation_kernel(
    propagation::config cfg,
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

    int gid = threadIdx.x + blockIdx.x * blockDim.x;
    /*if (gid >= navigator_states_view.size()) {
        return;
    }*/

    vecmem::device_vector<typename stepper_t::state> stepper_states(
        stepper_states_view);
    // vecmem::device_vector<unsigned int> stepper_results(stepper_res_view);
    vecmem::device_vector<typename navigator_t::state> navigation_states(
        navigator_states_view);
    // vecmem::device_vector<unsigned int>
    // navigation_results(navigator_res_view);

    if (gid < navigator_states_view.size()) {
        const auto dir = stepper_states.at(gid)().dir();
        printf("%f: [%f, %f, %f]\n", navigation_states.at(gid)(), dir[0],
               dir[1], dir[2]);
    }

    // bar[0] and bar[1] track if buffers buffer_0 and buffer_1 are ready to be
    // filled, while bar[2] and bar[3] track if buffers buffer_0 and buffer_1
    // are filled-in respectively
    __shared__ barrier_t bar[4];

    auto block = cooperative_groups::this_thread_block();
    if (block.thread_rank() < 4)
        init(bar + block.thread_rank(), block.size());
    block.sync();

    if (block.thread_rank() < 4 * warpSize) {
        printf("Thread %d doing stepping\n", gid);
        do_step(bar, bar + 2, navigation_states, stepper_states, gid);
    } else {
        printf("Thread %d doing navigation\n", gid);
        navigate(bar, bar + 2, stepper_states, navigation_states, gid);
    }

    // Continuous stepping
    /*if (block.thread_rank() < warpSize) {
        do_step<barrier_t>(stepper_states, stepper_results, navigation_states,
    navigation_results, gid);
    }
    // Continuous navigation
    else if (warpSize <=
             block.thread_rank() && block.thread_rank() < warpSize) {
        navigate<barrier_t>(stepper_states, stepper_results, navigation_states,
    navigation_results, gid);
    }
    // Continous calls to user functor
    else {

        // auto actor_state_refs =
        // actor_chain_t::setup_actor_states(actor_states);

        // Particle hypothesis
        // auto &ptc = p_state._stepping.particle_hypothesis();
        // p_state.set_particle(update_particle_hypothesis(ptc,
        // tracks.at(gid)));

        // run_actors(*device_actor_state_ptr);
    }*/
}

template <typename actor_chain_t>
typename actor_chain_t::state_tuple *setup_actor_states(
    typename actor_chain_t::state_tuple *input_actor_states) {

    // Copy the actor state blueprint to the device
    using actor_state_t = typename actor_chain_t::state_tuple;
    actor_state_t *device_actor_state_ptr{nullptr};

    cudaError_t success =
        cudaMalloc((void **)&device_actor_state_ptr, sizeof(actor_state_t));
    assert(success == cudaSuccess);

    success = cudaMemcpy(device_actor_state_ptr, input_actor_states,
                         sizeof(actor_state_t), cudaMemcpyHostToDevice);
    assert(success == cudaSuccess);

    return device_actor_state_ptr;
}

template <typename actor_chain_t>
void release_actor_states(
    typename actor_chain_t::state_tuple *device_actor_state_ptr) {
    [[maybe_unused]] cudaError_t success = cudaFree(device_actor_state_ptr);
    assert(success == cudaSuccess);
}

template <typename navigator_t, typename stepper_t, typename actor_chain_t>
void run_propagation_init_kernel(
    const propagation::config &cfg,
    typename navigator_t::detector_type::view_type det_view,
    typename stepper_t::magnetic_field_type field_view,
    vecmem::data::vector_view<typename stepper_t::free_track_parameters_type>
        tracks_view,
    vecmem::data::vector_view<typename stepper_t::state> stepper_states_view,
    vecmem::data::vector_view<typename navigator_t::state>
        navigator_states_view,
    vecmem::data::vector_view<typename actor_chain_t::state_tuple>
        actor_states_view,
    typename actor_chain_t::state_tuple *input_actor_states) {

    constexpr int thread_dim = 256;
    int block_dim = (tracks_view.size() + thread_dim - 1) / thread_dim;

    // Copy blueprint actor states to device
    auto *device_actor_state_ptr =
        setup_actor_states<actor_chain_t>(input_actor_states);

    // run the test kernel
    propagation_init_kernel<navigator_t, stepper_t, actor_chain_t>
        <<<block_dim, thread_dim>>>(cfg, det_view, field_view, tracks_view,
                                    stepper_states_view, navigator_states_view,
                                    actor_states_view, device_actor_state_ptr);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    // All states are set up: blueprint no longer needed
    release_actor_states<actor_chain_t>(device_actor_state_ptr);
}

template <typename navigator_t, typename stepper_t, typename actor_chain_t>
void run_propagation_kernel(
    const propagation::config &cfg,
    vecmem::data::vector_view<typename stepper_t::state> stepper_states_view,
    vecmem::data::vector_view<unsigned int> stepper_res_view,
    vecmem::data::vector_view<typename navigator_t::state>
        navigator_states_view,
    vecmem::data::vector_view<unsigned int> navigator_res_view,
    vecmem::data::vector_view<typename actor_chain_t::state_tuple>
        actor_states_view) {

    std::cout << "Number of tracks " << stepper_states_view.size() << std::endl;

    constexpr int thread_dim = 256;
    int block_dim = (stepper_states_view.size() + thread_dim - 1) / thread_dim;

    // run the propagation loop
    propagation_kernel<navigator_t, stepper_t, actor_chain_t>
        <<<block_dim, thread_dim>>>(cfg, stepper_states_view,
                                    navigator_states_view, actor_states_view);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/// Macro declaring the template instantiations for the different detector types
#define DECLARE_PROPAGATOR(METADATA, CHAIN, FIELD)                           \
                                                                             \
    template void run_propagation_init_kernel<                               \
        navigator_type<METADATA>, stepper_type<METADATA, FIELD>,             \
        CHAIN<detector<METADATA>::algebra_type>>(                            \
        const propagation::config &, typename detector<METADATA>::view_type, \
        covfie::field_view<FIELD>,                                           \
        vecmem::data::vector_view<                                           \
            free_track_parameters<detector<METADATA>::algebra_type>>,        \
        vecmem::data::vector_view<                                           \
            typename stepper_type<METADATA, FIELD>::state>,                  \
        vecmem::data::vector_view<typename navigator_type<METADATA>::state>, \
        vecmem::data::vector_view<                                           \
            typename CHAIN<detector<METADATA>::algebra_type>::state_tuple>,  \
        typename CHAIN<detector<METADATA>::algebra_type>::state_tuple *);    \
                                                                             \
    template void run_propagation_kernel<                                    \
        navigator_type<METADATA>, stepper_type<METADATA, FIELD>,             \
        CHAIN<detector<METADATA>::algebra_type>>(                            \
        const propagation::config &,                                         \
        vecmem::data::vector_view<                                           \
            typename stepper_type<METADATA, FIELD>::state>,                  \
        vecmem::data::vector_view<unsigned int>,                             \
        vecmem::data::vector_view<typename navigator_type<METADATA>::state>, \
        vecmem::data::vector_view<unsigned int>,                             \
        vecmem::data::vector_view<                                           \
            typename CHAIN<detector<METADATA>::algebra_type>::state_tuple>);

DECLARE_PROPAGATOR(benchmarks::default_metadata, empty_chain, const_field_t)
DECLARE_PROPAGATOR(benchmarks::default_metadata, default_chain, const_field_t)

DECLARE_PROPAGATOR(benchmarks::toy_metadata, empty_chain, const_field_t)
DECLARE_PROPAGATOR(benchmarks::toy_metadata, default_chain, const_field_t)

}  // namespace detray
