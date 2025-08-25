/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/propagation_config.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/log.hpp"

// Detray test include(s)
#include "detray/benchmarks/types.hpp"
#include "detray/test/common/bfield.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// System include(s)
#include <algorithm>
#include <cassert>
#include <chrono>
#include <ctime>
#include <iostream>
#include <ratio>
#include <string>

namespace detray {

// Define propagator type
template <concepts::algebra algebra_t>
using empty_chain = actor_chain<>;

template <concepts::algebra algebra_t>
using default_chain = actor_chain<parameter_transporter<algebra_t>,
                                  pointwise_material_interactor<algebra_t>,
                                  parameter_resetter<algebra_t>>;

using const_field_t = bfield::const_bknd_t<benchmarks::scalar>;

template <typename metadata_t, typename bfield_t>
using stepper_type = rk_stepper<covfie::field_view<bfield_t>,
                                typename detector<metadata_t>::algebra_type>;

template <typename metadata_t>
using navigator_type =
    caching_navigator<detector<metadata_t, device_container_types>>;

/// Propagation that state aggregates a stepping and a navigation state. It
/// also keeps references to the actor states.
template <typename navigator_t, typename stepper_t>
struct propagation_state {
    typename stepper_t::state _stepping;
    typename navigator_t::state _navigation;
    typename navigator_t::detector_type::geometry_context _context{};
};

/// Setup all of the data states on device
///
/// @param cfg the propagation configuration
/// @param det_view the detector vecmem view
/// @param field_data the magentic field view (maybe an empty field)
/// @param tracks_data the track collection view
/// @param navigation_cache_view the navigation cache vecemem view
/// @param opt which propagation to run (sync vs. unsync)
template <typename navigator_t, typename stepper_t, typename actor_chain_t>
void run_propagation_init_kernel(
    const propagation::config *, typename navigator_t::detector_type *,
    typename navigator_t::detector_type::view_type,
    typename stepper_t::magnetic_field_type,
    vecmem::data::vector_view<typename stepper_t::free_track_parameters_type>,
    vecmem::data::vector_view<propagation_state<navigator_t, stepper_t>>,
    vecmem::data::vector_view<typename actor_chain_t::state_tuple>,
    typename actor_chain_t::state_tuple *);

/// Launch the propagation kernel
///
/// @param cfg the propagation configuration
/// @param det_view the detector vecmem view
/// @param field_data the magentic field view (maybe an empty field)
/// @param tracks_data the track collection view
/// @param navigation_cache_view the navigation cache vecemem view
/// @param opt which propagation to run (sync vs. unsync)
template <typename navigator_t, typename stepper_t, typename actor_chain_t>
void run_propagation_kernel(
    const propagation::config *, const typename navigator_t::detector_type *,
    vecmem::data::vector_view<propagation_state<navigator_t, stepper_t>>,
    vecmem::data::vector_view<unsigned int>,
    vecmem::data::vector_view<unsigned int>,
    vecmem::data::vector_view<typename actor_chain_t::state_tuple>);

/// Copy the propagation config to device
typename propagation::config *setup_config(
    const typename propagation::config *);

/// Release the config allocation after the propagation is done
void release_config(typename propagation::config *);

/// Allocate space for the device detector opbject to be shared between kenels
template <typename device_detector_t>
device_detector_t *setup_device_detector(typename device_detector_t::view_type);

/// Release device detector
template <typename device_detector_t>
void release_device_detector(device_detector_t *);

/// Copy a blueprint actor state to device, which will be taken and set up
/// correctly for every track in the propagation init kernel
/// @returns the device
template <typename actor_chain_t>
typename actor_chain_t::state_tuple *setup_actor_states(
    typename actor_chain_t::state_tuple *);

/// Release the blueprint actor allocation after all actor states are set up
template <typename actor_chain_t>
void release_actor_states(typename actor_chain_t::state_tuple *);

/// Device Propagation benchmark
template <typename navigator_t, typename stepper_t, typename actor_chain_t>
struct cuda_propagation {
    /// Detector dependent types
    using metadata_t = typename navigator_t::detector_type::metadata;
    using detector_t = detector<metadata_t>;
    using device_detector_t = detector<metadata_t, device_container_types>;
    using bfield_view_t = typename stepper_t::magnetic_field_type;
    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    using state = propagation_state<navigator_t, stepper_t>;

    /// The propagation configuration
    propagation::config m_cfg{};

    /// Default construction
    cuda_propagation() = default;

    /// Construct from an externally provided configuration @param cfg
    explicit cuda_propagation(const propagation::config &cfg) : m_cfg{cfg} {}

    /// @return the propagation configuration
    propagation::config &config() { return m_cfg; }

    /// Prepare data and run propagation loop
    template <typename bfield_bkn_t>
    inline void operator()(
        vecmem::memory_resource *dev_mr, const detector_t *det,
        const bfield_bkn_t *field,
        dvector<free_track_parameters<algebra_t>> *tracks,
        typename actor_chain_t::state_tuple *input_actor_states) const {

        assert(dev_mr != nullptr);
        assert(tracks != nullptr);
        assert(det != nullptr);
        assert(field != nullptr);
        assert(input_actor_states != nullptr);

        // Helper object for performing memory copies (to CUDA devices)
        vecmem::cuda::copy cuda_cpy;

        // Copy the detector to device and get its view
        auto det_buffer = detray::get_buffer(*det, *dev_mr, cuda_cpy);
        auto det_view = detray::get_data(det_buffer);

        bfield_view_t field_view(*field);

        // Copy the track collection to device
        auto track_buffer =
            detray::get_buffer(vecmem::get_data(*tracks), *dev_mr, cuda_cpy);

        // Launch the propagator test for GPU device
        propagate(dev_mr, det_view, field_view, track_buffer,
                  input_actor_states);
    }

    /// Run propagation loop
    inline void operator()(
        vecmem::memory_resource *dev_mr, const detector_t::view_type det_view,
        const bfield_view_t field_view,
        vecmem::data::vector_view<free_track_parameters<algebra_t>> tracks_view,
        typename actor_chain_t::state_tuple *input_actor_states) const {

        // Launch the propagator test for GPU device
        propagate(dev_mr, det_view, field_view, tracks_view,
                  input_actor_states);
    }

    private:
    /// Prepare data and run propagation loop
    inline void propagate(
        vecmem::memory_resource *dev_mr, const detector_t::view_type det_view,
        const bfield_view_t field_view,
        vecmem::data::vector_view<free_track_parameters<algebra_t>> tracks_view,
        typename actor_chain_t::state_tuple *input_actor_states) const {

        assert(dev_mr != nullptr);
        assert(tracks_view.size() != 0u);
        assert(input_actor_states != nullptr);

        // Helper object for performing memory copies (to CUDA devices)
        vecmem::cuda::copy cuda_cpy;

        const unsigned int n_tracks{
            static_cast<unsigned int>(tracks_view.size())};

        // Buffer for the stepper states
        vecmem::data::vector_buffer<state> stepper_states_buffer(n_tracks,
                                                                 *dev_mr);
        cuda_cpy.setup(stepper_states_buffer)->wait();
        auto stepper_states_view = vecmem::get_data(stepper_states_buffer);

        // Buffer for the actor states
        vecmem::data::vector_buffer<typename actor_chain_t::state_tuple>
            actor_states_buffer(n_tracks, *dev_mr);
        cuda_cpy.setup(actor_states_buffer)->wait();
        auto actor_states_view = vecmem::get_data(actor_states_buffer);

        // Copy config to device
        auto *pinned_cfg_ptr = setup_config(&m_cfg);
        // Allocate memory for device detector object (not detector data!)
        auto *pinned_detector_ptr =
            setup_device_detector<device_detector_t>(det_view);

        // Copy blueprint actor states to device
        auto *pinned_actor_state_ptr =
            setup_actor_states<actor_chain_t>(input_actor_states);

        std::chrono::high_resolution_clock::time_point t1_init =
            std::chrono::high_resolution_clock::now();
        // Initialize the propagation (fill the state buffers)
        run_propagation_init_kernel<navigator_t, stepper_t, actor_chain_t>(
            pinned_cfg_ptr, pinned_detector_ptr, det_view, field_view,
            tracks_view, stepper_states_view, actor_states_view,
            pinned_actor_state_ptr);
        std::chrono::high_resolution_clock::time_point t2_init =
            std::chrono::high_resolution_clock::now();

        // All states are set up: blueprint no longer needed
        release_actor_states<actor_chain_t>(pinned_actor_state_ptr);

        const auto init_time =
            std::chrono::duration_cast<std::chrono::duration<double>>(t2_init -
                                                                      t1_init);
        const double init_ms{init_time.count() * 1000.};

        // Navigation results
        vecmem::data::vector_buffer<unsigned int> stepper_res_buffer(n_tracks,
                                                                     *dev_mr);
        cuda_cpy.setup(stepper_res_buffer)->wait();
        auto stepper_res_view = vecmem::get_data(stepper_res_buffer);

        // Navigation results
        vecmem::data::vector_buffer<unsigned int> navigation_res_buffer(
            n_tracks, *dev_mr);
        cuda_cpy.setup(navigation_res_buffer)->wait();
        auto navigation_res_view = vecmem::get_data(navigation_res_buffer);

        std::chrono::high_resolution_clock::time_point t1_prop =
            std::chrono::high_resolution_clock::now();
        // Launch the propagation for GPU device
        run_propagation_kernel<navigator_t, stepper_t, actor_chain_t>(
            pinned_cfg_ptr, pinned_detector_ptr, stepper_states_view,
            stepper_res_view, navigation_res_view, actor_states_view);
        std::chrono::high_resolution_clock::time_point t2_prop =
            std::chrono::high_resolution_clock::now();

        // Deallocate the propagation config and detector
        release_config(pinned_cfg_ptr);
        release_device_detector(pinned_detector_ptr);

        const auto prop_time =
            std::chrono::duration_cast<std::chrono::duration<double>>(t2_prop -
                                                                      t1_prop);
        const double prop_ms{prop_time.count() * 1000.};

        DETRAY_INFO_HOST("Propagator took: "
                         << init_ms + prop_ms << "ms ("
                         << (init_ms + prop_ms) / (n_tracks / 3000.)
                         << " ms/evt)");

        DETRAY_INFO_HOST("  -> Init: " << init_ms << "ms ("
                                       << init_ms / (n_tracks / 3000.)
                                       << " ms/evt)");

        DETRAY_INFO_HOST("  -> Propagation: " << prop_ms << "ms ("
                                              << prop_ms / (n_tracks / 3000.)
                                              << " ms/evt)\n");
    }
};

}  // namespace detray
