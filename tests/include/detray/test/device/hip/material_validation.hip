/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/detail/hip_definitions.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "material_validation.hpp"

namespace detray::hip {

template <typename detector_t>
__global__ void material_validation_kernel(
    typename detector_t::view_type det_data, const propagation::config cfg,
    vecmem::data::vector_view<
        free_track_parameters<typename detector_t::algebra_type>>
        tracks_view,
    vecmem::data::vector_view<
        material_validator::material_record<typename detector_t::scalar_type>>
        mat_records_view,
    vecmem::data::jagged_vector_view<
        material_validator::material_params<typename detector_t::scalar_type>>
        mat_steps_view) {

    using detector_device_t =
        detector<typename detector_t::metadata, device_container_types>;
    using algebra_t = typename detector_device_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;

    using stepper_t = line_stepper<algebra_t>;
    using navigator_t = navigator<detector_device_t>;
    // Propagator with full covariance transport, pathlimit aborter and
    // material tracer
    using material_tracer_t =
        material_validator::material_tracer<scalar_t, vecmem::device_vector>;
    using pathlimit_aborter_t = pathlimit_aborter<scalar_t>;
    using actor_chain_t =
        actor_chain<pathlimit_aborter_t, parameter_transporter<algebra_t>,
                    parameter_resetter<algebra_t>,
                    pointwise_material_interactor<algebra_t>,
                    material_tracer_t>;
    using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

    detector_device_t det(det_data);

    vecmem::device_vector<free_track_parameters<algebra_t>> tracks(tracks_view);
    vecmem::device_vector<typename material_tracer_t::material_record_type>
        mat_records(mat_records_view);
    vecmem::jagged_device_vector<
        typename material_tracer_t::material_params_type>
        mat_steps(mat_steps_view);

    int trk_id = threadIdx.x + blockIdx.x * blockDim.x;
    if (trk_id >= tracks.size()) {
        return;
    }

    propagator_t p{cfg};

    // Create the actor states
    typename pathlimit_aborter_t::state aborter_state{cfg.stepping.path_limit};
    typename pointwise_material_interactor<algebra_t>::state interactor_state{};
    typename material_tracer_t::state mat_tracer_state{mat_steps.at(trk_id)};

    auto actor_states =
        ::detray::tie(aborter_state, interactor_state, mat_tracer_state);

    // Run propagation
    typename navigator_t::state::view_type nav_view{};
    typename propagator_t::state propagation(tracks[trk_id], det, nav_view);

    p.propagate(propagation, actor_states);

    // Record the accumulated material
    assert(mat_records.size() == tracks.size());
    mat_records.at(trk_id) = mat_tracer_state.get_material_record();
}

/// Launch the device kernel
template <typename detector_t>
void material_validation_device(
    typename detector_t::view_type det_view, const propagation::config &cfg,
    vecmem::data::vector_view<
        free_track_parameters<typename detector_t::algebra_type>> &tracks_view,
    vecmem::data::vector_view<
        material_validator::material_record<typename detector_t::scalar_type>>
        &mat_records_view,
    vecmem::data::jagged_vector_view<
        material_validator::material_params<typename detector_t::scalar_type>>
        &mat_steps_view) {

    constexpr int thread_dim = 2 * WARP_SIZE;
    int block_dim = tracks_view.size() / thread_dim + 1;

    // run the test kernel
    hipLaunchKernelGGL((material_validation_kernel<detector_t>) , dim3(block_dim) , dim3(thread_dim) , 0 , 0 ,
        det_view, cfg, tracks_view, mat_records_view, mat_steps_view);


    // hip error check
    DETRAY_HIP_ERROR_CHECK(hipGetLastError());
    DETRAY_HIP_ERROR_CHECK(hipDeviceSynchronize());
}

/// Macro declaring the template instantiations for the different detector types
#define DECLARE_MATERIAL_VALIDATION(METADATA)                                 \
                                                                              \
    template void material_validation_device<detector<METADATA>>(             \
        typename detector<METADATA>::view_type, const propagation::config &,  \
        vecmem::data::vector_view<                                            \
            free_track_parameters<typename detector<METADATA>::algebra_type>> \
            &,                                                                \
        vecmem::data::vector_view<material_validator::material_record<        \
            typename detector<METADATA>::scalar_type>> &,                     \
        vecmem::data::jagged_vector_view<material_validator::material_params< \
            typename detector<METADATA>::scalar_type>> &);

DECLARE_MATERIAL_VALIDATION(test::default_metadata)
DECLARE_MATERIAL_VALIDATION(test::toy_metadata)
DECLARE_MATERIAL_VALIDATION(test::default_telescope_metadata)

}  // namespace detray::hip
