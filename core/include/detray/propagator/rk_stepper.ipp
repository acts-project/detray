/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/geometry/tracking_volume.hpp"

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE inline void
detray::rk_stepper<magnetic_field_t, algebra_t, constraint_t, policy_t,
                   inspector_t, array_t>::state::advance_track() {

    const auto& sd = this->_step_data;
    const scalar_type h{this->_step_size};
    const scalar_type h_6{h * static_cast<scalar_type>(1. / 6.)};
    auto& track = this->_track;
    auto pos = track.pos();
    auto dir = track.dir();

    // Update the track parameters according to the equations of motion
    // Reference: Eq (82) of https://doi.org/10.1016/0029-554X(81)90063-X
    pos = pos + h * (sd.t[0u] + h_6 * (sd.dtds[0] + sd.dtds[1] + sd.dtds[2]));
    track.set_pos(pos);

    // Reference: Eq (82) of https://doi.org/10.1016/0029-554X(81)90063-X
    dir =
        dir + h_6 * (sd.dtds[0] + 2.f * (sd.dtds[1] + sd.dtds[2]) + sd.dtds[3]);
    dir = vector::normalize(dir);
    track.set_dir(dir);

    auto qop = track.qop();
    if (!(this->_mat == nullptr)) {
        // Reference: Eq (82) of https://doi.org/10.1016/0029-554X(81)90063-X
        qop =
            qop + h_6 * (sd.dqopds[0u] + 2.f * (sd.dqopds[1u] + sd.dqopds[2u]) +
                         sd.dqopds[3u]);
    }
    track.set_qop(qop);

    // Update path length
    this->_path_length += h;
    this->_abs_path_length += math::fabs(h);
    this->_s += h;
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE inline void detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::advance_jacobian(const detray::stepping::config& cfg) {
    /// The calculations are based on ATL-SOFT-PUB-2009-002. The update of the
    /// Jacobian matrix is requires only the calculation of eq. 17 and 18.
    /// Since the terms of eq. 18 are currently 0, this matrix is not needed
    /// in the calculation. The matrix A from eq. 17 consists out of 3
    /// different parts. The first one is given by the upper left 3x3 matrix
    /// that are calculated by the derivatives dF/dT (called dFdT) and dG/dT
    /// (calles dGdT). The second is given by the top 3 lines of the rightmost
    /// column. This is calculated by dFdqop and dGdqop. The remaining non-zero
    /// term is calculated directly. The naming of the variables is explained in
    /// eq. 11 and are directly related to the initial problem in eq. 7. The
    /// evaluation is based by propagating the parameters T and lambda as given
    /// in eq. 16 and evaluating the derivations for matrix A.
    /// @note The translation for u_{n+1} in eq. 7 is in this case a
    /// 3-dimensional vector without a dependency of Lambda or lambda neither in
    /// u_n nor in u_n'. The second and fourth eq. in eq. 14 have the constant
    /// offset matrices h * Id and Id respectively. This involves that the
    /// constant offset does not exist for rectangular matrix dGdu' (due to the
    /// missing Lambda part) and only exists for dFdu' in dlambda/dlambda.

    // Set transport matrix (D) and update Jacobian transport
    //( JacTransport = D * JacTransport )
    auto D = matrix_operator().template identity<e_free_size, e_free_size>();

    const auto& sd = this->_step_data;
    const scalar_type h{this->_step_size};
    // const auto& mass = this->_mass;
    auto& track = this->_track;

    // Half step length
    const scalar_type h2{h * h};
    const scalar_type half_h{h * 0.5f};
    const scalar_type h_6{h * static_cast<scalar_type>(1. / 6.)};

    // 3X3 Identity matrix
    const matrix_type<3, 3> I33 = matrix_operator().template identity<3, 3>();

    // Initialize derivatives
    std::array<matrix_type<3u, 3u>, 4u> dkndt{I33, I33, I33, I33};
    std::array<vector3_type, 4u> dkndqop;
    std::array<matrix_type<3u, 3u>, 4u> dkndr;
    std::array<scalar_type, 4u> dqopn_dqop{1.f, 1.f, 1.f, 1.f};

    /*---------------------------------------------------------------------------
     *  dk_n/dt1
     *    = qop_n * (dt_n/dt1 X B_n)
     *      + qop_n * ( t_n X dB_n/dt1 ),
     *  where dB_n/dt1 == dB_n/dr_n * dr_n/dt1.
     *
     *  The second term is non-zero only for inhomogeneous magnetic fields
     *
     *  Note that [ t_n = t1 + h * d(t_{n-1})/ds) ] as indicated by
     *  Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X

     *  [ Table for dt_n/dt1 ]
     *  dt1/dt1 = I
     *  dt2/dt1 = d( t1 + h/2 * dt1/ds ) / dt1 = I + h/2 * dk1/dt1
     *  dt3/dt1 = d( t1 + h/2 * dt2/ds ) / dt1 = I + h/2 * dk2/dt1
     *  dt4/dt1 = d( t1 + h * dt3/ds ) / dt1 = I + h * dk3/dt1
     *
     *  [ Table for dr_n/dt1 ]
     *  dr1/dt1 = 0
     *  dr2/dt1 = d(r1 + h/2 * t1 + h^2/8 dt1/ds)/dt1 = h/2 * I + h^2/8 dk1/dt1
     *  dr3/dt1 = d(r1 + h/2 * t1 + h^2/8 dt1/ds)/dt1 = h/2 * I + h^2/8 dk1/dt1
     *  dr4/dt1 = d(r1 + h * t1 + h^2/2 dt3/ds)/dt1 = h * I + h^2/2 dk3/dt1
     *
     *  Note that
     *  d/dr [ F(T) X B ]  = dF(T)/dr (X) B, where (X) means the column wise
     *  cross product
    ---------------------------------------------------------------------------*/

    /*---------------------------------------------------------------------------
     *  dk_n/dqop_1
     *    = dqop_n/dqop1 * ( t_n X B_n )
     *      + qop_n * ( dt_n/dqop1 X B_n )
     *      + qop_n * ( t_n X dB_n/dqop1 ),
     *  where dB_n/dqop1 = dB_n/dr_n * dr_n/dqop1
     *
     *  Note that [ qop_n = qop1 + h * dqop_{n-1}/ds) ] as indicated by
     *  Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
     *
     *  [ Table for dqop_n/dqop1 ]
     *  dqop1/dqop1 = 1
     *  dqop2/dqop1 = 1 + h/2 * d(dqop1/ds)/dqop1
     *  dqop3/dqop1 = 1 + h/2 * d(dqop2/ds)/dqop1
     *  dqop4/dqop1 = 1 + h * d(dqop3/ds)/dqop1
     *
     *  [ Table for dt_n/dqop1 ]
     *  dt1/dqop1 = 0
     *  dt2/dqop1 = d(t1 + h/2 dt1/ds)/dqop1 = h/2 * dk1/dqop1
     *  dt3/dqop1 = d(t1 + h/2 dt2/ds)/dqop1 = h/2 * dk2/dqop1
     *  dt4/dqop1 = d(t1 + h dt3/ds)/dqop1 = h * dk3/dqop1
     *
     *  [ Table for dr_n/dqop1 ]
     *  dr1/dqop1 = 0
     *  dr2/dqop1 = d(r1 + h/2 * t1 + h^2/8 dt1/ds)/dqop1 = h^2/8 * dk1/dqop1
     *  dr3/dqop1 = d(r1 + h/2 * t1 + h^2/8 dt1/ds)/dqop1 = h^2/8 * dk1/dqop1
     *  dr4/dqop1 = d(r1 + h * t1 + h^2/2 dt3/ds)/dqop1 = h^2/2 dk3/dqop1
    ---------------------------------------------------------------------------*/

    /*---------------------------------------------------------------------------
     *  dk_n/dr1
     *    = qop_n * ( dt_n/dr1 X B_n )
     *      + qop_n * ( t_n X dB_n/dr1 ),
     *  where dB_n/dr1 = dB_n/dr_n * dr_n/dr1
     *
     *  [ Table for dt_n/dr1 ]
     *  dt1/dr1 = 0
     *  dt2/dr1 = d(t1 + h/2 * dt1/ds)/dr1 = h/2 * dk1/dr1
     *  dt2/dr1 = d(t1 + h/2 * dt2/ds)/dr1 = h/2 * dk2/dr1
     *  dt3/dr1 = d(t1 + h * dt3/ds)/dr1 = h * dk3/dr1
     *
     *  [ Table for dr_n/dr1 ]
     *  dr1/dr1 = I
     *  dr2/dr1 = (r1 + h/2 * t1 + h^2/8 dt1/ds ) / dr1 = I + h^2/8 dk1/dr1
     *  dr3/dr1 = (r1 + h/2 * t1 + h^2/8 dt1/ds ) / dr1 = I + h^2/8 dk1/dr1
     *  dr4/dr1 = (r1 + h * t1 + h^2/2 dt3/ds ) / dr1 = I + h^2/2 dk3/dr1
    ---------------------------------------------------------------------------*/

    /*---------------------------------------------------------------------------
     *  d(dqop_n/ds)/dqop1
     *
     *  Useful equation:
     *  dqop/ds = qop^3 * E * (-dE/ds) / q^2 = - qop^3 * E g / q^2

     *  [ Table for d(dqop_n/ds)/dqop1 ]
     *  d(dqop1/ds)/dqop1 = dqop1/ds * (1/qop * (3 - p^2/E^2) + 1/g1 * dg1dqop1
     *  d(dqop2/ds)/dqop1 = d(dqop2/ds)/dqop2 * (1 + h/2 * d(dqop1/ds)/dqop1)
     *  d(dqop3/ds)/dqop1 = d(dqop3/ds)/dqop3 * (1 + h/2 * d(dqop2/ds)/dqop1)
     *  d(dqop4/ds)/dqop1 = d(dqop4/ds)/dqop4 * (1 + h * d(dqop3/ds)/dqop1)
    ---------------------------------------------------------------------------*/

    if (!cfg.use_eloss_gradient) {
        getter::element(D, e_free_qoverp, e_free_qoverp) = 1.f;
    } else {
        // Pre-calculate dqop_n/dqop1
        const scalar_type d2qop1dsdqop1 = this->d2qopdsdqop(sd.qop[0u]);

        dqopn_dqop[0u] = 1.f;
        dqopn_dqop[1u] = 1.f + half_h * d2qop1dsdqop1;

        const scalar_type d2qop2dsdqop1 =
            this->d2qopdsdqop(sd.qop[1u]) * dqopn_dqop[1u];
        dqopn_dqop[2u] = 1.f + half_h * d2qop2dsdqop1;

        const scalar_type d2qop3dsdqop1 =
            this->d2qopdsdqop(sd.qop[2u]) * dqopn_dqop[2u];
        dqopn_dqop[3u] = 1.f + h * d2qop3dsdqop1;

        const scalar_type d2qop4dsdqop1 =
            this->d2qopdsdqop(sd.qop[3u]) * dqopn_dqop[3u];

        /*-----------------------------------------------------------------
         * Calculate the first terms of d(dqop_n/ds)/dqop1
        -------------------------------------------------------------------*/

        getter::element(D, e_free_qoverp, e_free_qoverp) =
            1.f + h_6 * (d2qop1dsdqop1 + 2.f * (d2qop2dsdqop1 + d2qop3dsdqop1) +
                         d2qop4dsdqop1);
    }

    /*-----------------------------------------------------------------
     * Calculate the first terms of dk_n/dt1
    -------------------------------------------------------------------*/
    // dk1/dt1
    dkndt[0u] =
        sd.qop[0u] * mat_helper().column_wise_cross(dkndt[0u], sd.b_first);

    // dk2/dt1
    dkndt[1u] = dkndt[1u] + half_h * dkndt[0u];
    dkndt[1u] =
        sd.qop[1u] * mat_helper().column_wise_cross(dkndt[1u], sd.b_middle);

    // dk3/dt1
    dkndt[2u] = dkndt[2u] + half_h * dkndt[1u];
    dkndt[2u] =
        sd.qop[2u] * mat_helper().column_wise_cross(dkndt[2u], sd.b_middle);

    // dk4/dt1
    dkndt[3u] = dkndt[3u] + h * dkndt[2u];
    dkndt[3u] =
        sd.qop[3u] * mat_helper().column_wise_cross(dkndt[3u], sd.b_last);

    /*-----------------------------------------------------------------
     * Calculate the first and second terms of dk_n/dqop1
    -------------------------------------------------------------------*/
    // dk1/dqop1
    dkndqop[0u] = dqopn_dqop[0u] * vector::cross(sd.t[0u], sd.b_first);

    // dk2/dqop1
    dkndqop[1u] = dqopn_dqop[1u] * vector::cross(sd.t[1u], sd.b_middle) +
                  sd.qop[1u] * half_h * vector::cross(dkndqop[0u], sd.b_middle);

    // dk3/dqop1
    dkndqop[2u] = dqopn_dqop[2u] * vector::cross(sd.t[2u], sd.b_middle) +
                  sd.qop[2u] * half_h * vector::cross(dkndqop[1u], sd.b_middle);

    // dk4/dqop1
    dkndqop[3u] = dqopn_dqop[3u] * vector::cross(sd.t[3u], sd.b_last) +
                  sd.qop[3u] * h * vector::cross(dkndqop[2u], sd.b_last);

    // Calculate dkndr in case of considering B field gradient
    if (cfg.use_field_gradient) {

        // Positions and field gradients at initial, middle and final points of
        // the fourth order RKN
        vector3_type r_ini = track.pos();
        vector3_type r_mid =
            r_ini + half_h * sd.t[0u] + h2 * 0.125f * sd.dtds[0u];
        vector3_type r_fin = r_ini + h * sd.t[0u] + h2 * 0.5f * sd.dtds[2u];

        matrix_type<3, 3> dBdr_ini = evaluate_field_gradient(r_ini);
        matrix_type<3, 3> dBdr_mid = evaluate_field_gradient(r_mid);
        matrix_type<3, 3> dBdr_fin = evaluate_field_gradient(r_fin);

        /*-----------------------------------------------------------------
         * Calculate all terms of dk_n/dr1
        -------------------------------------------------------------------*/

        // dk1/dr1
        dkndr[0u] =
            -sd.qop[0u] * mat_helper().column_wise_cross(dBdr_ini, sd.t[0u]);

        // dk2/dr1
        dkndr[1u] = sd.qop[1u] * mat_helper().column_wise_cross(
                                     half_h * dkndr[0u], sd.b_middle);
        dkndr[1u] = dkndr[1u] -
                    sd.qop[1u] * mat_helper().column_wise_cross(
                                     dBdr_mid * (I33 + h2 * 0.125 * dkndr[0u]),
                                     sd.t[1u]);

        // dk3/dr1
        dkndr[2u] = sd.qop[2u] * mat_helper().column_wise_cross(
                                     half_h * dkndr[1u], sd.b_middle);
        dkndr[2u] = dkndr[2u] -
                    sd.qop[2u] * mat_helper().column_wise_cross(
                                     dBdr_mid * (I33 + h2 * 0.125 * dkndr[0u]),
                                     sd.t[2u]);

        // dk4/dr1
        dkndr[3u] = sd.qop[3u] *
                    mat_helper().column_wise_cross(h * dkndr[2u], sd.b_last);
        dkndr[3u] =
            dkndr[3u] -
            sd.qop[3u] * mat_helper().column_wise_cross(
                             dBdr_fin * (I33 + h2 * 0.5 * dkndr[2u]), sd.t[3u]);

        // Set dF/dr1 and dG/dr1
        auto dFdr = matrix_operator().template identity<3, 3>();
        auto dGdr = matrix_operator().template identity<3, 3>();
        dFdr = dFdr + h * h_6 * (dkndr[0u] + dkndr[1u] + dkndr[2u]);
        dGdr = h_6 * (dkndr[0u] + 2.f * (dkndr[1u] + dkndr[2u]) + dkndr[3u]);

        matrix_operator().set_block(D, dFdr, 0u, 0u);
        matrix_operator().set_block(D, dGdr, 4u, 0u);
    }

    // Set dF/dt1 and dG/dt1
    auto dFdt = matrix_operator().template identity<3, 3>();
    auto dGdt = matrix_operator().template identity<3, 3>();
    dFdt = dFdt + h_6 * (dkndt[0u] + dkndt[1u] + dkndt[2u]);
    dFdt = h * dFdt;
    dGdt = dGdt + h_6 * (dkndt[0u] + 2.f * (dkndt[1u] + dkndt[2u]) + dkndt[3u]);

    matrix_operator().set_block(D, dFdt, 0u, 4u);
    matrix_operator().set_block(D, dGdt, 4u, 4u);

    // Set dF/dqop1 and dG/dqop1
    vector3_type dFdqop = h * h_6 * (dkndqop[0u] + dkndqop[1u] + dkndqop[2u]);
    vector3_type dGdqop =
        h_6 * (dkndqop[0u] + 2.f * (dkndqop[1u] + dkndqop[2u]) + dkndqop[3u]);
    matrix_operator().set_block(D, dFdqop, 0u, 7u);
    matrix_operator().set_block(D, dGdqop, 4u, 7u);

    this->_jac_transport = D * this->_jac_transport;
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE inline auto detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::evaluate_dqopds(const std::size_t i, const scalar_type h,
                                     const scalar_type dqopds_prev,
                                     const detray::stepping::config& cfg)
    -> scalar_type {

    const auto& track = this->_track;
    const scalar_type qop = track.qop();
    auto& sd = this->_step_data;

    if (this->_mat == nullptr) {
        sd.qop[i] = qop;
        return 0.f;
    } else {

        if (cfg.use_mean_loss) {
            if (i == 0u) {
                sd.qop[i] = qop;
            } else {

                // qop_n is calculated recursively like the direction of
                // evaluate_dtds.
                //
                // https://doi.org/10.1016/0029-554X(81)90063-X says:
                // "For y  we  have  similar  formulae  as  for x, for y' and
                // \lambda similar  formulae as for  x'"
                sd.qop[i] = qop + h * dqopds_prev;
            }
        }
        return this->dqopds(sd.qop[i]);
    }
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE inline auto detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::evaluate_dtds(const vector3_type& b_field,
                                   const std::size_t i, const scalar_type h,
                                   const vector3_type& dtds_prev,
                                   const scalar_type qop) -> vector3_type {
    auto& track = this->_track;
    const auto dir = track.dir();
    auto& sd = this->_step_data;

    if (i == 0u) {
        sd.t[i] = dir;
    } else {
        // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
        sd.t[i] = dir + h * dtds_prev;
    }

    // dtds = qop * (t X B) from Lorentz force
    return qop * vector::cross(sd.t[i], b_field);
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE inline auto detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::evaluate_field_gradient(const point3_type& pos)
    -> matrix_type<3, 3> {

    matrix_type<3, 3> dBdr = matrix_operator().template zero<3, 3>();

    constexpr auto delta{1e-1f * unit<scalar_type>::mm};

    for (unsigned int i = 0; i < 3; i++) {

        point3_type dpos1 = pos;
        dpos1[i] += delta;
        const auto bvec1_tmp =
            this->_magnetic_field.at(dpos1[0], dpos1[1], dpos1[2]);
        vector3_type bvec1;
        bvec1[0u] = bvec1_tmp[0u];
        bvec1[1u] = bvec1_tmp[1u];
        bvec1[2u] = bvec1_tmp[2u];

        point3_type dpos2 = pos;
        dpos2[i] -= delta;
        const auto bvec2_tmp =
            this->_magnetic_field.at(dpos2[0], dpos2[1], dpos2[2]);
        vector3_type bvec2;
        bvec2[0u] = bvec2_tmp[0u];
        bvec2[1u] = bvec2_tmp[1u];
        bvec2[2u] = bvec2_tmp[2u];

        const vector3_type gradient = (bvec1 - bvec2) * (1.f / (2.f * delta));

        getter::element(dBdr, 0u, i) = gradient[0u];
        getter::element(dBdr, 1u, i) = gradient[1u];
        getter::element(dBdr, 2u, i) = gradient[2u];
    }

    return dBdr;
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE inline auto
detray::rk_stepper<magnetic_field_t, algebra_t, constraint_t, policy_t,
                   inspector_t, array_t>::state::dtds() const -> vector3_type {

    // In case there was no step before
    if (this->_path_length == 0.f) {
        const point3_type pos = this->_track.pos();

        const auto bvec_tmp = this->_magnetic_field.at(pos[0], pos[1], pos[2]);
        vector3_type bvec;
        bvec[0u] = bvec_tmp[0u];
        bvec[1u] = bvec_tmp[1u];
        bvec[2u] = bvec_tmp[2u];

        return this->_track.qop() * vector::cross(this->_track.dir(), bvec);
    }
    return this->_step_data.dtds[3u];
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE inline auto
detray::rk_stepper<magnetic_field_t, algebra_t, constraint_t, policy_t,
                   inspector_t, array_t>::state::dqopds() const -> scalar_type {

    // In case there was no step before
    if (this->_path_length == 0.f) {
        return this->dqopds(this->_track.qop());
    }

    return this->_step_data.dqopds[3u];
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::dqopds(const scalar_type qop) const -> scalar_type {

    // d(qop)ds is zero for empty space
    if (this->_mat == nullptr) {
        return 0.f;
    }

    const auto& mat = this->volume_material();

    const scalar_type q = this->_ptc.charge();
    const scalar_type p = q / qop;
    const scalar_type mass = this->_ptc.mass();
    const scalar_type E = math::sqrt(p * p + mass * mass);

    // Compute stopping power
    const scalar_type stopping_power =
        interaction<scalar_type>().compute_stopping_power(mat, this->_ptc,
                                                          {mass, qop, q});

    // Assert that a momentum is a positive value
    assert(p >= 0.f);

    // d(qop)ds, which is equal to (qop) * E * (-dE/ds) / p^2
    // or equal to (qop)^3 * E * (-dE/ds) / q^2
    return qop * qop * qop * E * stopping_power / (q * q);
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::d2qopdsdqop(const scalar_type qop) const -> scalar_type {

    if (this->_mat == nullptr) {
        return 0.f;
    }

    const auto& mat = this->volume_material();
    const scalar_type q = this->_ptc.charge();
    const scalar_type p = q / qop;
    const scalar_type p2 = p * p;

    const auto& mass = this->_ptc.mass();
    const scalar_type E2 = p2 + mass * mass;

    // Interaction object
    interaction<scalar_type> I;

    // g = dE/ds = -1 * (-dE/ds) = -1 * stopping power
    const detail::relativistic_quantities<scalar_type> rq(mass, qop, q);
    const scalar_type g = -1.f * I.compute_stopping_power(mat, this->_ptc, rq);

    // dg/d(qop) = -1 * derivation of stopping power
    const scalar_type dgdqop =
        -1.f * I.derive_stopping_power(mat, this->_ptc, rq);

    // d(qop)/ds = - qop^3 * E * g / q^2
    const scalar_type dqopds = this->dqopds(qop);

    // Check Eq 3.12 of
    // (https://iopscience.iop.org/article/10.1088/1748-0221/4/04/P04016/meta)
    return dqopds * (1.f / qop * (3.f - p2 / E2) + 1.f / g * dgdqop);
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
template <typename propagation_state_t>
DETRAY_HOST_DEVICE inline bool detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::step(propagation_state_t& propagation,
                   const detray::stepping::config& cfg) const {

    // Get stepper and navigator states
    state& stepping = propagation._stepping;
    auto& magnetic_field = stepping._magnetic_field;
    auto& navigation = propagation._navigation;

    if (stepping._step_size == 0.f) {
        stepping._step_size = cfg.min_stepsize;
    } else if (stepping._step_size > 0) {
        stepping._step_size = math::min(stepping._step_size, navigation());
    } else {
        stepping._step_size = math::max(stepping._step_size, navigation());
    }

    const point3_type pos = stepping().pos();

    auto vol = navigation.get_volume();
    if (vol.has_material()) {
        stepping._mat = vol.material_parameters(pos);
    } else {
        stepping._mat = nullptr;
    }

    auto& sd = stepping._step_data;

    scalar_type error_estimate{0.f};

    // First Runge-Kutta point
    const auto bvec = magnetic_field.at(pos[0], pos[1], pos[2]);
    sd.b_first[0] = bvec[0];
    sd.b_first[1] = bvec[1];
    sd.b_first[2] = bvec[2];

    // qop should be recalcuated at every point
    // Reference: Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
    sd.dqopds[0u] = stepping.evaluate_dqopds(0u, 0.f, 0.f, cfg);
    sd.dtds[0u] = stepping.evaluate_dtds(
        sd.b_first, 0u, 0.f, vector3_type{0.f, 0.f, 0.f}, sd.qop[0u]);

    const auto estimate_error = [&](const scalar_type& h) -> scalar {
        // State the square and half of the step size
        const scalar_type h2{h * h};
        const scalar_type half_h{h * 0.5f};

        // Second Runge-Kutta point
        // qop should be recalcuated at every point
        // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
        const point3_type pos1 =
            pos + half_h * sd.t[0u] + h2 * 0.125f * sd.dtds[0u];
        const auto bvec1 = magnetic_field.at(pos1[0], pos1[1], pos1[2]);
        sd.b_middle[0] = bvec1[0];
        sd.b_middle[1] = bvec1[1];
        sd.b_middle[2] = bvec1[2];

        sd.dqopds[1u] =
            stepping.evaluate_dqopds(1u, half_h, sd.dqopds[0u], cfg);
        sd.dtds[1u] = stepping.evaluate_dtds(sd.b_middle, 1u, half_h,
                                             sd.dtds[0u], sd.qop[1u]);

        // Third Runge-Kutta point
        // qop should be recalcuated at every point
        // Reference: Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
        sd.dqopds[2u] =
            stepping.evaluate_dqopds(2u, half_h, sd.dqopds[1u], cfg);
        sd.dtds[2u] = stepping.evaluate_dtds(sd.b_middle, 2u, half_h,
                                             sd.dtds[1u], sd.qop[2u]);

        // Last Runge-Kutta point
        // qop should be recalcuated at every point
        // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
        const point3_type pos2 = pos + h * sd.t[0u] + h2 * 0.5f * sd.dtds[2u];
        const auto bvec2 = magnetic_field.at(pos2[0], pos2[1], pos2[2]);
        sd.b_last[0] = bvec2[0];
        sd.b_last[1] = bvec2[1];
        sd.b_last[2] = bvec2[2];

        sd.dqopds[3u] = stepping.evaluate_dqopds(3u, h, sd.dqopds[2u], cfg);
        sd.dtds[3u] =
            stepping.evaluate_dtds(sd.b_last, 3u, h, sd.dtds[2u], sd.qop[3u]);

        // Compute and check the local integration error estimate
        // @Todo
        constexpr const scalar_type one_sixth{
            static_cast<scalar_type>(1. / 6.)};
        const vector3_type err_vec =
            one_sixth * h2 *
            (sd.dtds[0u] - sd.dtds[1u] - sd.dtds[2u] + sd.dtds[3u]);
        error_estimate = getter::norm(err_vec);

        return error_estimate;
    };

    scalar_type error{1e20f};

    // Whenever navigator::init() is called the step size is set to navigation
    // path length (navigation()). We need to reduce it down to make error small
    // enough
    for (unsigned int i_t = 0u; i_t < cfg.max_rk_updates; i_t++) {
        stepping.count_trials();

        error = math::max(estimate_error(stepping._step_size),
                          static_cast<scalar_type>(1e-20));

        // Error is small enough
        // ---> break and advance track
        if (error <= 4.f * cfg.rk_error_tol) {
            break;
        }
        // Error estimate is too big
        // ---> Make step size smaller and esimate error again
        else {

            scalar_type step_size_scaling =
                math::sqrt(math::sqrt(cfg.rk_error_tol / error));

            stepping._step_size *= step_size_scaling;

            // Run inspection while the stepsize is getting adjusted
            stepping.run_inspector(cfg, "Adjust stepsize: ", i_t + 1,
                                   step_size_scaling);
        }
    }

    // Update navigation direction
    const step::direction step_dir = stepping._step_size >= 0.f
                                         ? step::direction::e_forward
                                         : step::direction::e_backward;
    stepping.set_direction(step_dir);

    // Check constraints
    if (math::fabs(stepping._step_size) >
        math::fabs(
            stepping.constraints().template size<>(stepping.direction()))) {

        // Run inspection before step size is cut
        stepping.run_inspector(cfg, "Before constraint: ");

        stepping.set_step_size(
            stepping.constraints().template size<>(stepping.direction()));
    }

    // Advance track state
    stepping.advance_track();

    // Advance jacobian transport
    if (cfg.do_covariance_transport) {
        stepping.advance_jacobian(cfg);
    }

    // Call navigation update policy
    typename rk_stepper::policy_type{}(stepping.policy_state(), propagation);

    const auto step_size_scaling = static_cast<scalar_type>(
        math::min(math::max(math::sqrt(math::sqrt(cfg.rk_error_tol / error)),
                            static_cast<scalar_type>(0.25)),
                  static_cast<scalar_type>(4.)));

    // Save the current step size
    stepping._prev_step_size = stepping._step_size;

    // Update the step size
    stepping._step_size *= step_size_scaling;

    // Run final inspection
    stepping.run_inspector(cfg, "Step complete: ");

    return true;
}
