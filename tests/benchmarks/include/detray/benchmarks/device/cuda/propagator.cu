/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/benchmarks/device/cuda/propagator.hpp"
#include "detray/definitions/detail/cuda_definitions.hpp"
#include "detray/utils/log.hpp"

// CUDA include(s)
#include <cooperative_groups.h>

#include <cuda/barrier>

namespace detray {

template <typename navigator_t, typename stepper_t, typename actor_chain_t>
struct propagation_queue {
    using barrier_type = ::cuda::barrier<::cuda::thread_scope_block>;
    using data_buffer_type =
        std::array<propagation_state<navigator_t, stepper_t>, 2u>;
    using slot_type = std::array<int, 2>;

    barrier_type *m_barriers;

    // Thread group that works on this queue
    cooperative_groups::thread_group &m_group;

    // Data to work on, shared between all blocks!
    data_buffer_type& m_data;

    DETRAY_DEVICE propagation_queue(
        data_buffer_type& data,
        cooperative_groups::thread_group &thread_group, barrier_type barriers[], const int bar_idx)
        : m_data{data}, m_group{thread_group}, m_barriers{barriers} {

        if (m_group.thread_rank() < 4) {
            init(m_barriers + bar_idx + m_group.thread_rank(), m_group.size());
        }
        m_group.sync();
    }

    propagation_queue(const propagation_queue &data) = default;

    /// Set the thread group that will work on this queue
    DETRAY_DEVICE void thread_group(cooperative_groups::thread_group group) {
        m_group = group;
    }

    /// @returns the thread group that works on this queue
    DETRAY_DEVICE cooperative_groups::thread_group thread_group() const {
        return m_group;
    }

    DETRAY_DEVICE constexpr std::size_t size() const { return m_data.size(); }

    DETRAY_DEVICE constexpr const auto &data() const { return m_data; }

    DETRAY_DEVICE constexpr int slot_index(const int step) const {
        return step % 2;
    }

    DETRAY_DEVICE constexpr auto& get_data(const int step) const {
        return m_data[slot_index(step)];
    }

    DETRAY_DEVICE constexpr auto& get_data(const int step) {
        return m_data[slot_index(step)];
    }

    DETRAY_DEVICE constexpr decltype(auto) get_stepping_barrier(
        const int step) {
        return m_barriers[slot_index(step)];
    }

    DETRAY_DEVICE constexpr decltype(auto) get_navigation_barrier(
        const int step) {
        return m_barriers[slot_index(step) + 2];
    }

    DETRAY_DEVICE constexpr const auto &navigation_states(
        const int step) const {
        return get_data(step)._navigation;
    }

    DETRAY_DEVICE constexpr auto &navigation_states(const int step) {
        return get_data(step)._navigation;
    }

    DETRAY_DEVICE constexpr const auto &stepper_states(const int step) const {
        return get_data(step)._stepping;
    }

    DETRAY_DEVICE constexpr auto &stepper_states(const int step) {
        return get_data(step)._stepping;
    }

    DETRAY_DEVICE constexpr bool is_alive(const int step) const {
        return navigation_states(step).is_alive();
    }

    DETRAY_DEVICE inline void wait_for_navigation(const int step) {
        DETRAY_INFO_HOST_DEVICE(
            "Thread %ld (Step): Wait for navigation (slot: %d, iteration %d)",
            m_group.thread_rank(), slot_index(step), step);
        // Wait until the navigation triggers stepping
        __syncwarp();
        get_stepping_barrier(step).arrive_and_wait();
    }

    DETRAY_DEVICE inline void wait_for_stepping(const int step) {
        DETRAY_INFO_HOST_DEVICE(
            "Thread %ld (Nav): Wait for stepping (slot: %d, iteration %d)",
            m_group.thread_rank(), slot_index(step), step);
        // Wait until the stepper triggers navigation
        __syncwarp();
        get_navigation_barrier(step).arrive_and_wait();
    }

    DETRAY_DEVICE inline barrier_type::arrival_token trigger_navigation(
        const int step) {
        DETRAY_INFO_HOST_DEVICE(
            "Thread %ld (Step): Trigger navigation (slot: %d, iteration %d)",
            m_group.thread_rank(), slot_index(step), step);
        __syncwarp();
        return get_navigation_barrier(step).arrive();
    }

    DETRAY_DEVICE inline barrier_type::arrival_token trigger_stepping(
        const int step) {
        DETRAY_INFO_HOST_DEVICE(
            "Thread %ld (Nav): Trigger stepper (slot: %d, iteration %d)",
            m_group.thread_rank(), slot_index(step), step);
        __syncwarp();
        return get_stepping_barrier(step).arrive();
    }
};

template <typename collection_t>
__device__ std::size_t scan(collection_t &c) {

    for (std::size_t i = 0; i < c.size(); ++i) {
    }
}

/// Run the stepper
template <typename navigator_t, typename stepper_t, typename actor_chain_t>
__device__ inline void take_step(
    const stepping::config *cfg,
    propagation_queue<navigator_t, stepper_t, actor_chain_t> &queue,
    [[maybe_unused]] int thread_idx) {

    using barrier_t = ::cuda::barrier<::cuda::thread_scope_block>;

    using detector_t = typename navigator_t::detector_type;
    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;

    using stepper_state_t = typename stepper_t::state;
    using navigation_state_t = typename navigator_t::state;
    using propagation_state_t = propagation_state<navigator_t, stepper_t>;

    constexpr stepper_t stepper{};
    constexpr bool reset_stepsize{true};

    assert(cfg != nullptr);

    for (int i = 0; i < 10; ++i) {

        queue.wait_for_navigation(i);

        if (queue.is_alive(i)) {
            const navigation_state_t &navigation = queue.navigation_states(i);
            stepper_state_t &stepping = queue.stepper_states(i);

            DETRAY_INFO_HOST_DEVICE(
                "Thread %d (Step): Take step ready (iteration %d)", thread_idx,
                i);

            const auto &track = stepping();

            const auto vol = navigation.current_volume();
            const material<scalar_t> *vol_mat_ptr =
                vol.has_material() ? vol.material_parameters(track.pos())
                                   : nullptr;

            stepper.step(navigation(), stepping, *cfg, reset_stepsize,
                         vol_mat_ptr);

            typename stepper_t::policy_type{}(stepping.policy_state(),
                                              queue.get_data(i));
            DETRAY_INFO_HOST_DEVICE(
                "Thread %d (Step): Step length %f (iteration %d)", thread_idx,
                stepping.step_size(), i);
        }

        barrier_t::arrival_token token = queue.trigger_navigation(i);
    }

    DETRAY_INFO_HOST_DEVICE("Thread %d, block %d (Step): Stepping finished",
                        thread_idx, blockIdx.x);
}

/// Run the navigator
template <typename navigator_t, typename stepper_t, typename actor_chain_t>
__device__ inline void navigate(
    const navigation::config *cfg,
    propagation_queue<navigator_t, stepper_t, actor_chain_t> &queue,
    [[maybe_unused]] int thread_idx) {

    using barrier_t = ::cuda::barrier<::cuda::thread_scope_block>;

    using detector_t = typename navigator_t::detector_type;
    using stepper_state_t = typename stepper_t::state;
    using navigation_state_t = typename navigator_t::state;
    using propagation_state_t = propagation_state<navigator_t, stepper_t>;

    assert(cfg != nullptr);
    // assert(static_cast<int>(queue.size()) > 512 * blockIdx.x);

    // Block-local thread index (adjusted for the start of the warp
    // specialization)
    // const int thread_idx =
    //    cooperative_groups::this_thread_block().thread_rank() - blockDim.x /
    //    2;

    barrier_t::arrival_token token1 = queue.trigger_stepping(0);
    barrier_t::arrival_token token2 = queue.trigger_stepping(1);

    for (int i = 0; i < 10; ++i) {

        queue.wait_for_stepping(i);

        if (queue.is_alive(i)) {
            navigation_state_t &navigation = queue.navigation_states(i);
            const stepper_state_t &stepper_state = queue.stepper_states(i);

            DETRAY_INFO_HOST_DEVICE(
                "Thread %d (Nav): Run navigation (iteration % d) ", thread_idx,
                i);

            constexpr navigator_t navigator{};
            typename detector_t::geometry_context gctx{};
            navigator.update(stepper_state(), navigation, *cfg, gctx);
        }

        barrier_t::arrival_token token = queue.trigger_stepping(i);
    }

    DETRAY_INFO_HOST_DEVICE("Thread %d, block %d (Nav): Navigation finished",
                        thread_idx, blockIdx.x);
}

/// Specialize the warps to run stepping, navigation and actors independently
template <typename navigator_t, typename stepper_t, typename actor_chain_t>
__global__ void __launch_bounds__(256, 4) propagation_kernel(
    const propagation::config *pinned_cfg,
    const typename navigator_t::detector_type *pinned_detector_ptr,
    typename stepper_t::magnetic_field_type field_view,
    vecmem::data::vector_view<typename stepper_t::free_track_parameters_type>
        tracks_view,
    // vecmem::data::vector_view<unsigned int> stepper_res_view,
    // vecmem::data::vector_view<unsigned int> navigator_res_view,
    typename actor_chain_t::state_tuple *pinned_actor_state_ptr) {

    using detector_device_t = typename navigator_t::detector_type;
    using algebra_t = typename detector_device_t::algebra_type;

    using track_t = typename stepper_t::free_track_parameters_type;
    using stepper_state_t = typename stepper_t::state;
    using navigation_state_t = typename navigator_t::state;
    using propagation_state_t = propagation_state<navigator_t, stepper_t>;

    const vecmem::device_vector<track_t> tracks(tracks_view);

    constexpr navigator_t navigator{};
    constexpr actor_chain_t run_actors{};

    typename detector_device_t::geometry_context gctx{};

    // Setup some small amount of data in shared memory
    __shared__ stepping::config step_cfg[1];
    __shared__ navigation::config nav_cfg[1];
    __shared__ typename actor_chain_t::state_tuple actor_states_sh[2];

    auto block = cooperative_groups::this_thread_block();
    if (block.thread_rank() == 0) {
        step_cfg[0] = pinned_cfg->stepping;
        nav_cfg[0] = pinned_cfg->navigation;
        // Create the actor states on a fresh copy
        actor_states_sh[0] = *pinned_actor_state_ptr;
        actor_states_sh[1] = *pinned_actor_state_ptr;
    }
    block.sync();

    /// Every thread initializes one set of states
    assert(blockDim.y == blockDim.z == 1);

    int gid = threadIdx.x + blockIdx.x * blockDim.x;
    if (2*gid + 1 >= tracks.size()) {
        return;
    }

    ///
    /// Run propagation init
    ///

    // Propagation states: Two per track
    std::array<propagation_state_t, 2u> propagation{propagation_state_t{stepper_state_t{tracks.at(2*gid), field_view}, navigation_state_t{*pinned_detector_ptr}, gctx}, propagation_state_t{stepper_state_t{tracks.at(2*gid + 1), field_view}, navigation_state_t{*pinned_detector_ptr}, gctx}};

    for (int i = 0; i < 2; ++i) {
        stepper_state_t& stepping = propagation[i]._stepping;
        navigation_state_t& navigation = propagation[i]._navigation;

        // Initialize the navigation
        navigator.init(stepping(), navigation, *nav_cfg, gctx);
        assert(navigation.is_alive());

        // Run all registered actors/aborters
        auto actor_states_ref =
            actor_chain_t::setup_actor_states(actor_states_sh[i]);
        run_actors(actor_states_ref, propagation[i]);
        assert(!stepping().is_invalid());
        assert(!stepping.bound_params().is_invalid());

        // Update the navigation information, in case the actors changed the track
        navigator.update(stepping(), navigation, *nav_cfg, gctx);

        assert(navigation.is_alive());
        assert(std::isfinite(navigation()));

        DETRAY_INFO_HOST_DEVICE("Thread %d (Init): Dist to next %f mm",
                            block.thread_rank(),
                            navigation());
    }

    ///
    /// Run propagation steps
    ///
    using barrier_t = ::cuda::barrier<::cuda::thread_scope_block>;

    // bar[0] and bar[1] track if buffers buffer_0 and buffer_1 are ready to be
    // filled, while bar[2] and bar[3] track if buffers buffer_0 and buffer_1
    // are filled-in respectively
    __shared__ barrier_t barriers[4 * 256 / 64];

    cooperative_groups::thread_group warp = cooperative_groups::tiled_partition<32>(cooperative_groups::this_thread_block());

    if (block.thread_rank() < blockDim.x / 2) {
        DETRAY_INFO_HOST_DEVICE("Thread %d (%d) in block %d doing stepping", gid,
                            block.thread_rank(), blockIdx.x);

        const int bar_idx{4 * (block.thread_rank() % 32)};
        assert(bar_idx < 4 * 256 / 64);

        propagation_queue<navigator_t, stepper_t, actor_chain_t> queue{
            propagation, warp, barriers, bar_idx};

        take_step(step_cfg, queue, gid);
    } else {
        DETRAY_INFO_HOST_DEVICE("Thread %d (%d) in block %d doing navigation", gid,
                            block.thread_rank(), blockIdx.x);

        const int bar_idx{4 * ((block.thread_rank() - blockDim.x / 2 ) % 32)};
        assert(bar_idx < 4 * 256 / 64);

        propagation_queue<navigator_t, stepper_t, actor_chain_t> queue{
            propagation, warp, barriers, bar_idx};

        navigate(nav_cfg, queue, gid);
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
                                       
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    return pinned_config_ptr;
}

void release_config(typename propagation::config *pinned_config_ptr) {
    DETRAY_CUDA_ERROR_CHECK(cudaFreeHost(pinned_config_ptr));
}

template <typename device_detector_t>
device_detector_t *setup_device_detector(
    typename device_detector_t::view_type det_view) {

    // Build a device detector type (the interal pointers and capacities refer
    // to the already allocated vecmem device buffers)
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
    const propagation::config *pinned_cfg,
    const typename navigator_t::detector_type *pinned_detector_ptr,
    typename stepper_t::magnetic_field_type field_view,
    vecmem::data::vector_view<typename stepper_t::free_track_parameters_type>
        tracks_view,
    vecmem::data::vector_view<unsigned int> stepper_res_view,
    vecmem::data::vector_view<unsigned int> navigator_res_view,
    typename actor_chain_t::state_tuple *pinned_actor_state_ptr) {

    int thread_dim = math::min(256, static_cast<int>(tracks_view.size()));
    // One block handles 256 tracks (2 per thread)
    int block_dim = static_cast<int>(
        (tracks_view.size() / 2 + thread_dim - 1) / thread_dim);

    DETRAY_INFO_HOST_DEVICE("# Tracks: %ld", tracks_view.size());
    DETRAY_INFO_HOST_DEVICE("# threads per block: %d", thread_dim);
    DETRAY_INFO_HOST_DEVICE("# blocks: %d", block_dim);
    DETRAY_INFO_HOST_DEVICE("# threads: %d", thread_dim * block_dim);

    // run the propagation loop
    propagation_kernel<navigator_t, stepper_t, actor_chain_t>
        <<<block_dim, thread_dim>>>(pinned_cfg, pinned_detector_ptr, field_view, tracks_view, pinned_actor_state_ptr);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/// Macro declaring the template instantiations for the different detector types
#define DECLARE_PROPAGATOR(METADATA, CHAIN, FIELD)                          \
                                                                            \
    template void run_propagation_kernel<                                   \
        navigator_type<METADATA>, stepper_type<METADATA, FIELD>,            \
        CHAIN<detector<METADATA>::algebra_type>>(                           \
        const propagation::config *,                                        \
        const detector<METADATA, device_container_types> *,                 \
        covfie::field_view<FIELD>,                                          \
        vecmem::data::vector_view<                                          \
            free_track_parameters<detector<METADATA>::algebra_type>>,       \
        vecmem::data::vector_view<unsigned int>,                            \
        vecmem::data::vector_view<unsigned int>,                            \
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
