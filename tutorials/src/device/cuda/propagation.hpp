// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Project include(s).
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"

// Tutorial include(s)
#include "detray/tutorial/types.hpp"

// Covfie include(s)
#include <covfie/cuda/backend/primitive/cuda_device_array.hpp>

namespace detray::tutorial {

// Detector
using metadata_t = detray::tutorial::toy_metadata;
using detector_host_t = detector<metadata_t, host_container_types>;
using detector_device_t = detector<metadata_t, device_container_types>;

using algebra_t = metadata_t::algebra_type;
using scalar = detray::tutorial::scalar;

namespace bfield::cuda {

// Inhomogeneous field (cuda)
using inhom_bknd_t = covfie::backend::affine<covfie::backend::linear<
    covfie::backend::strided<covfie::vector::vector_d<std::size_t, 3>,
                             covfie::backend::cuda_device_array<
                                 covfie::vector::vector_d<scalar, 3>>>>>;

}  // namespace bfield::cuda

// Navigator
using navigator_t = navigator<detector_device_t>;
using intersection_t = navigator_t::intersection_type;

// Stepper
using host_field_t = covfie::field<detray::bfield::inhom_bknd_t<scalar>>;
using device_field_t =
    covfie::field<detray::tutorial::bfield::cuda::inhom_bknd_t>;
using stepper_t = rk_stepper<device_field_t::view_t, algebra_t>;

// Actors
using actor_chain_t =
    actor_chain<pathlimit_aborter<scalar>, parameter_transporter<algebra_t>,
                pointwise_material_interactor<algebra_t>,
                parameter_resetter<algebra_t>>;

// Propagator
using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

/// Propagation tutorial function
void propagation(
    typename detector_host_t::view_type det_data,
    typename device_field_t::view_t field_data,
    const vecmem::data::vector_view<free_track_parameters<algebra_t>>
        tracks_data);

}  // namespace detray::tutorial
