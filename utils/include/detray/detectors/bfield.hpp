/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/io/common/read_bfield.hpp"

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/nearest_neighbour.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/field.hpp>

namespace detray::bfield {

/// Constant bfield (host and device)
using const_bknd_t =
    covfie::backend::constant<covfie::vector::vector_d<scalar, 3>,
                              covfie::vector::vector_d<scalar, 3>>;

using const_field_t = covfie::field<const_bknd_t>;

/// Inhomogeneous field (host)
using inhom_bknd_t = covfie::backend::affine<
    covfie::backend::nearest_neighbour<covfie::backend::strided<
        covfie::vector::ulong3,
        covfie::backend::array<covfie::vector::vector_d<scalar, 3>>>>>;

using inhom_field_t = covfie::field<inhom_bknd_t>;

/// @returns a constant covfie field constructed from the field vector @param B
inline const_field_t create_const_field(
    const __plugin::vector3<detray::scalar> &B) {
    return const_field_t{covfie::make_parameter_pack(
        const_bknd_t::configuration_t{B[0], B[1], B[2]})};
}

/// @returns a constant covfie field constructed from the field vector @param B
inline inhom_field_t create_inhom_field() {
    return io::read_bfield<inhom_field_t>(
        !std::getenv("DETRAY_BFIELD_FILE") ? ""
                                           : std::getenv("DETRAY_BFIELD_FILE"));
}

}  // namespace detray::bfield
