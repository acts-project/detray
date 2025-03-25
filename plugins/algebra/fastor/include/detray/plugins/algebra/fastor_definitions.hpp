/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

 #pragma once

 // Algebra-Plugins include
 #include "algebra/fastor_fastor.hpp"
 
 namespace detray {
 
 /// The plugin definition
 template <algebra::concepts::scalar scalar_t>
 using fastor = algebra::plugin::fastor<scalar_t>;
 
 namespace getter {
 
 using algebra::fastor::storage::block;
 using algebra::fastor::storage::element;
 using algebra::fastor::storage::set_block;
 using algebra::fastor::storage::vector;
 
 }  // namespace getter
 
 namespace vector {
 
 using algebra::fastor::math::cross;
 using algebra::fastor::math::dot;
 using algebra::fastor::math::eta;
 using algebra::fastor::math::norm;
 using algebra::fastor::math::normalize;
 
 using algebra::fastor::math::perp;
 using algebra::fastor::math::phi;
 using algebra::fastor::math::theta;
 
 }  // namespace vector
 
 namespace matrix {
 
 using algebra::fastor::math::determinant;
 using algebra::fastor::math::identity;
 using algebra::fastor::math::inverse;
 using algebra::fastor::math::set_identity;
 using algebra::fastor::math::set_zero;
 using algebra::fastor::math::transpose;
 using algebra::fastor::math::zero;
 
 }  // namespace matrix
 
 }  // namespace detray
 
