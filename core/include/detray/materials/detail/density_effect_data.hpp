/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/units.hpp"

// System include(s)
#include <climits>

namespace detray::detail {

template <typename scalar_t>
struct density_effect_data {

    using scalar_type = scalar_t;

    // Plasma energy in eV
    scalar_type m_plasma_energy = std::numeric_limits<scalar>::epsilon();
    // Sternheimer adjustment factor for the atomic excitation energies
    scalar_type m_rho = std::numeric_limits<scalar>::epsilon();
    // -C
    scalar_type m_nC = std::numeric_limits<scalar>::epsilon();
    // parameters in fitting formulas
    scalar_type m_X0 = std::numeric_limits<scalar>::epsilon();
    scalar_type m_X1 = std::numeric_limits<scalar>::epsilon();
    scalar_type m_m = std::numeric_limits<scalar>::epsilon();
    scalar_type m_a = std::numeric_limits<scalar>::epsilon();
    // Density-effect value delta(X_0)
    scalar_type m_delta0 = std::numeric_limits<scalar>::epsilon();
    // Upper bound for the error inherent in the fitting procedure
    scalar_type m_delta_max = std::numeric_limits<scalar>::epsilon();
    // Mean excitation energy in eV
    scalar_type m_I = std::numeric_limits<scalar>::epsilon();
};

// Null density
constexpr const density_effect_data<scalar> null_density{};

// H2 Liquid
constexpr const density_effect_data<scalar> H2_liquid_density{
    7.031, 1.546, 3.2632, 0.4759, 1.9215, 0.13483, 5.6249, 0., 0.021, 21.8};

// H2 Gas
constexpr const density_effect_data<scalar> H2_gas_density{
    0.263, 1.412, 9.5835, 1.8639, 3.2718, 0.14092, 5.7273, 0.0, 0.024, 19.2};

// Helium Gas
constexpr const density_effect_data<scalar> He_gas_density{
    0.263, 1.7, 11.1393, 2.2017, 3.6122, 0.13443, 5.8347, 0, 0.024, 41.8};

// Berylium
constexpr const density_effect_data<scalar> Be_density{
    26.096, 1.908, 2.7847, 0.0392, 1.6922, 0.80392, 2.4339, 0.14, 0.029, 63.7};

// N2 Gas
constexpr const density_effect_data<scalar> N2_gas_density{
    0.695, 1.984, 10.5400, 1.7378, 4.1323, 0.15349, 3.2125, 0.0, 0.086, 82.};

// O2 Gas
constexpr const density_effect_data<scalar> O2_gas_density{
    0.744, 2.314, 10.7004, 1.7541, 4.3213, 0.11778, 3.2913, 0.0, 0.101, 95.};

// Aluminium
constexpr const density_effect_data<scalar> Al_density{
    32.86, 2.18, 4.2395, 0.1708, 3.0127, 0.08024, 3.6345, 0.12, 0.061, 166.};

// Silicon
constexpr const density_effect_data<scalar> Si_density{
    31.055, 2.103, 4.4351, 0.2014, 2.8715, 0.14921, 3.2546, 0.14, 0.059, 173.};

// Argon Gas
constexpr const density_effect_data<scalar> Ar_gas_density{
    0.789, 1.753, 11.9480, 1.7635, 4.4855, 0.19714, 2.9618, 0.0, 0.037, 188.};

// Tungsten
constexpr const density_effect_data<scalar> W_density{
    80.315, 1.997, 5.4059, 0.2167, 3.496, 0.15509, 2.8447, 0.14, 0.027, 727.};

// Gold
constexpr const density_effect_data<scalar> Au_density{
    80.215, 1.926, 5.5747, 0.2021, 3.6979, 0.09756, 3.1101, 0.14, 0.020, 790.};

}  // namespace detray::detail