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

    density_effect_data() = default;

    constexpr density_effect_data(const scalar_type a, const scalar_type m,
                                  const scalar_type X0, const scalar_type X1,
                                  const scalar_type I, const scalar_type nC,
                                  const scalar_type delta0)
        : m_a(a),
          m_m(m),
          m_X0(X0),
          m_X1(X1),
          m_I(I * unit_constants::eV),
          m_nC(nC),
          m_delta0(delta0) {}

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE bool operator==(
        const density_effect_data<scalar_t> &rhs) const {
        return (m_a == rhs.get_A_density() && m_m == rhs.get_M_density() &&
                m_X0 == rhs.get_X0_density() && m_X1 == rhs.get_X1_density() &&
                m_I == rhs.get_mean_excitation_energy() &&
                m_nC == rhs.get_C_density() &&
                m_delta0 == rhs.get_delta0_density());
    }

    DETRAY_HOST_DEVICE
    constexpr scalar_type get_A_density() const { return m_a; }

    DETRAY_HOST_DEVICE
    constexpr scalar_type get_M_density() const { return m_m; }

    DETRAY_HOST_DEVICE
    constexpr scalar_type get_X0_density() const { return m_X0; }

    DETRAY_HOST_DEVICE
    constexpr scalar_type get_X1_density() const { return m_X1; }

    DETRAY_HOST_DEVICE
    constexpr scalar_type get_mean_excitation_energy() const { return m_I; }

    DETRAY_HOST_DEVICE
    constexpr scalar_type get_C_density() const { return m_nC; }

    DETRAY_HOST_DEVICE
    constexpr scalar_type get_delta0_density() const { return m_delta0; }

    // Fitting parameters of Eq. 33.7 of RPP 2018
    scalar_type m_a = std::numeric_limits<scalar>::epsilon();
    scalar_type m_m = std::numeric_limits<scalar>::epsilon();
    scalar_type m_X0 = std::numeric_limits<scalar>::epsilon();
    scalar_type m_X1 = std::numeric_limits<scalar>::epsilon();
    // Mean excitation energy in eV
    scalar_type m_I = std::numeric_limits<scalar>::epsilon();
    // -C
    scalar_type m_nC = std::numeric_limits<scalar>::epsilon();
    // Density-effect value delta(X_0)
    scalar_type m_delta0 = std::numeric_limits<scalar>::epsilon();
};

// Values taken from https://pdg.lbl.gov/2022/AtomicNuclearProperties for Muon
// dEdX and range

// Null density
constexpr const density_effect_data<scalar> null_density{};

// H2 Liquid
constexpr const density_effect_data<scalar> H2_liquid_density{
    0.1348, 5.6249, 0.4400, 1.8856, 21.8, 3.0977, 0.00};

// H2 Gas
constexpr const density_effect_data<scalar> H2_gas_density{
    0.1409, 5.7273, 1.8639, 3.2718, 19.2, 9.5834, 0.00};

// Helium Gas
constexpr const density_effect_data<scalar> He_gas_density{
    0.1344, 5.8347, 2.2017, 3.6122, 41.8, 11.1393, 0.00};

// Berylium
constexpr const density_effect_data<scalar> Be_density{
    0.8039, 2.4339, 0.0592, 1.6922, 63.7, 2.7847, 0.14};

// N2 Gas
constexpr const density_effect_data<scalar> N2_gas_density{
    0.1535, 3.2125, 1.7378, 4.1323, 82.0, 10.5400, 0.00};

// O2 Gas
constexpr const density_effect_data<scalar> O2_gas_density{
    0.1178, 3.2913, 1.7541, 4.3213, 95.0, 10.7004, 0.00};

// Aluminium
constexpr const density_effect_data<scalar> Al_density{
    0.0802, 3.6345, 0.1708, 3.0127, 166.0, 4.2395, 0.12};

// Silicon
constexpr const density_effect_data<scalar> Si_density{
    0.1492, 3.2546, 0.2015, 2.8716, 173.0, 4.4355, 0.14};

// Argon Gas
constexpr const density_effect_data<scalar> Ar_gas_density{
    0.1971, 2.9618, 1.7635, 4.4855, 188.0, 11.9480, 0.00};

// Tungsten
constexpr const density_effect_data<scalar> W_density{
    0.1551, 2.8447, 0.2167, 3.4960, 727.0, 5.4059, 0.14};

// Gold
constexpr const density_effect_data<scalar> Au_density{
    0.0976, 3.1101, 0.2021, 3.6979, 790.0, 5.5747, 0.14};

}  // namespace detray::detail