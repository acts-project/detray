/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/io/covfie/read_bfield.hpp"

// Covfie include(s)
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/linear.hpp>
#include <covfie/core/backend/transformer/nearest_neighbour.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/vector.hpp>

namespace detray::bfield {

/// Constant bfield (host and device)
template <typename T>
using const_bknd_t = covfie::backend::constant<covfie::vector::vector_d<T, 3>,
                                               covfie::vector::vector_d<T, 3>>;

template <typename T>
using const_field_t = covfie::field<const_bknd_t<T>>;

/// Inhomogeneous field (host)
template <typename T>
using inhom_bknd_t =
    covfie::backend::affine<covfie::backend::linear<covfie::backend::strided<
        covfie::vector::vector_d<std::size_t, 3>,
        covfie::backend::array<covfie::vector::vector_d<T, 3>>>>>;

/// Inhomogeneous field with nearest neighbor (host)
template <typename T>
using inhom_bknd_nn_t = covfie::backend::affine<
    covfie::backend::nearest_neighbour<covfie::backend::strided<
        covfie::vector::ulong3,
        covfie::backend::array<covfie::vector::vector_d<T, 3>>>>>;

template <typename T>
using inhom_field_t = covfie::field<inhom_bknd_t<T>>;

/// @returns a constant covfie field constructed from the field vector @param B
template <typename T, concepts::vector3D vector3_t>
inline const_field_t<T> create_const_field(const vector3_t &B) {
    return const_field_t<T>{covfie::make_parameter_pack(
        typename const_bknd_t<T>::configuration_t{B[0], B[1], B[2]})};
}

/// @returns a constant covfie field constructed from the field vector @param B
template <typename T>
inline inhom_field_t<T> create_inhom_field() {
    return io::read_bfield<inhom_field_t<T>>(
        !std::getenv("DETRAY_BFIELD_FILE") ? ""
                                           : std::getenv("DETRAY_BFIELD_FILE"));
}

}  // namespace detray::bfield
