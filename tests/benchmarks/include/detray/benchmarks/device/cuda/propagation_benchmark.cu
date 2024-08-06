/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/benchmarks/device/cuda/propagation_benchmark.hpp"
#include "detray/core/detector_metadata.hpp"
#include "detray/definitions/detail/cuda_definitions.hpp"
#include "detray/detectors/toy_metadata.hpp"

namespace detray {

template <typename propagator_t>
__global__ void __launch_bounds__(256, 4) propagator_benchmark_kernel(
    propagation::config cfg,
    typename propagator_t::detector_type::view_type det_view,
    typename propagator_t::stepper_type::magnetic_field_type field_view,
    vecmem::data::vector_view<
        free_track_parameters<typename propagator_t::algebra_type>>
        tracks_view,
    vecmem::data::jagged_vector_view<
        typename propagator_t::navigator_type::intersection_type>
        nav_cache_view,
    const propagate_option opt) {

    using detector_device_t =
        detector<typename propagator_t::detector_type::metadata,
                 device_container_types>;
    using algebra_t = typename detector_device_t::algebra_type;
    using propagator_device_t =
        propagator<typename propagator_t::stepper_type,
                   navigator<detector_device_t>,
                   typename propagator_t::actor_chain_type>;
    using intersection_t = typename propagator_device_t::intersection_type;

    detector_device_t det(det_view);
    vecmem::device_vector<free_track_parameters<algebra_t>> tracks(tracks_view);
    vecmem::jagged_device_vector<intersection_t> nav_cache(nav_cache_view);

    int gid = threadIdx.x + blockIdx.x * blockDim.x;
    if (gid >= nav_cache.size()) {
        return;
    }

    // Create propagator
    propagator_device_t p{cfg};

    typename parameter_transporter<algebra_t>::state transporter_state{};
    typename pointwise_material_interactor<algebra_t>::state interactor_state{};
    typename parameter_resetter<algebra_t>::state resetter_state{};

    // Create the actor states
    auto actor_states =
        tie(transporter_state, interactor_state, resetter_state);

    // Create the propagator state
    typename propagator_device_t::state p_state(tracks.at(gid), field_view, det,
                                                nav_cache.at(gid));

    // Run propagation
    if (opt == propagate_option::e_unsync) {
        p.propagate(p_state, actor_states);
    } else if (opt == propagate_option::e_sync) {
        p.propagate_sync(p_state, actor_states);
    }
}

template <typename propagator_t>
void run_propagation_kernel(
    const propagation::config& cfg,
    typename propagator_t::detector_type::view_type det_view,
    typename propagator_t::stepper_type::magnetic_field_type field_view,
    vecmem::data::vector_view<
        free_track_parameters<typename propagator_t::algebra_type>>
        tracks_view,
    vecmem::data::jagged_vector_view<
        typename propagator_t::navigator_type::intersection_type>&
        candidates_data,
    const propagate_option opt) {

    constexpr int thread_dim = 256;
    int block_dim =
        static_cast<int>(candidates_data.size() + thread_dim - 1) / thread_dim;

    // run the test kernel
    propagator_benchmark_kernel<propagator_t><<<block_dim, thread_dim>>>(
        cfg, det_view, field_view, tracks_view, candidates_data, opt);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/// Macro declaring the template instantiations for the different detector types
#define DECLARE_PROPAGATION_BENCHMARK(METADATA, FIELD)                  \
                                                                        \
    template void                                                       \
    run_propagation_kernel<cuda_propagator_type<METADATA, FIELD>>(      \
        const propagation::config&, detector<METADATA>::view_type,      \
        covfie::field_view<FIELD>,                                      \
        vecmem::data::vector_view<                                      \
            free_track_parameters<detector<METADATA>::algebra_type>>,   \
        vecmem::data::jagged_vector_view<                               \
            cuda_propagator_type<METADATA, FIELD>::intersection_type>&, \
        const propagate_option);

DECLARE_PROPAGATION_BENCHMARK(default_metadata, bfield::const_bknd_t)
DECLARE_PROPAGATION_BENCHMARK(toy_metadata, bfield::const_bknd_t)

}  // namespace detray
