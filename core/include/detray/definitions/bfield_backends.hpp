/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/nearest_neighbour.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>

namespace detray::bfield {

// Constant bfield (host and device)
using const_bknd_t =
    covfie::backend::constant<covfie::vector::vector_d<scalar, 3>,
                              covfie::vector::vector_d<scalar, 3>>;

// Inhomogeneous field (host)
using inhom_bknd_t = covfie::backend::affine<
    covfie::backend::nearest_neighbour<covfie::backend::strided<
        covfie::vector::ulong3,
        covfie::backend::array<covfie::vector::vector_d<scalar, 3>>>>>;

}  // namespace detray::bfield
