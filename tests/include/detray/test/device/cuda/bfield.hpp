
/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"

// Covfie include(s)
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/clamp.hpp>
#include <covfie/core/backend/transformer/linear.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/vector.hpp>
#include <covfie/cuda/backend/primitive/cuda_device_array.hpp>

namespace detray::bfield::cuda {

// Inhomogeneous field (cuda)
template <typename T>
using inhom_bknd_t = covfie::backend::affine<
    covfie::backend::linear<covfie::backend::clamp<covfie::backend::strided<
        covfie::vector::vector_d<std::size_t, 3>,
        covfie::backend::cuda_device_array<covfie::vector::vector_d<T, 3>>>>>>;

// Test that the type is a valid backend for a CUDA field
static_assert(covfie::concepts::field_backend<inhom_bknd_t<float>>,
              "inhom_bknd_t is not a valid CUDA field backend type");

}  // namespace detray::bfield::cuda
