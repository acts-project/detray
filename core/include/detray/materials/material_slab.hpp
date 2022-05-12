/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

template <typename material_t>
struct material_slab {
    using scalar_type = typename material_t::scalar_type;

    material(const material_t& material, scalar_type thickness)
        : _material(material),
          _thickness(thickness),
          _thickness_in_X0(thickness / material.X0()),
          _thickness_in_L0(thickness / material.L0()) {}

    /// Access the (average) material parameters.
    constexpr const material_t& material() const { return _material; }
    /// Return the thickness.
    constexpr scalar_type thickness() const { return _thickness; }
    /// Return the radiation length fraction.
    constexpr scalar_type thickness_in_X0() const { return _thickness_in_X0; }
    /// Return the nuclear interaction length fraction.
    constexpr scalar_type thickness_in_L0() const { return _thickness_in_L0; }

    private:
    material_t _material;
    scalar_type _thickness = 0.0f;
    scalar_type _thickness_in_X0 = 0.0f;
    scalar_type _thickness_in_L0 = 0.0f;
};

}  // namespace detray