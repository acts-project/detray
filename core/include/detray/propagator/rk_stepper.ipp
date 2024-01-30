/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/geometry/detector_volume.hpp"

// System include(s)
#include <cmath>

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
void detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                        inspector_t, array_t>::state::advance_track() {

    using scalar_t = typename transform3_t::scalar_type;

    const auto& sd = this->_step_data;
    const scalar_t h{this->_step_size};
    const scalar_t h_6{h * static_cast<scalar_t>(1. / 6.)};
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
    if (!(this->_mat == vacuum<scalar_t>())) {
        // Reference: Eq (82) of https://doi.org/10.1016/0029-554X(81)90063-X
        qop =
            qop + h_6 * (sd.dqopds[0u] + 2.f * (sd.dqopds[1u] + sd.dqopds[2u]) +
                         sd.dqopds[3u]);
    }
    track.set_qop(qop);

    // Update path length
    this->_path_length += h;
    this->_s += h;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
void detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, inspector_t,
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

    using scalar_t = typename transform3_t::scalar_type;

    const auto& sd = this->_step_data;
    const scalar_t h{this->_step_size};
    // const auto& mass = this->_mass;
    auto& track = this->_track;

    // Half step length
    const scalar_t h2{h * h};
    const scalar_t half_h{h * 0.5f};
    const scalar_t h_6{h * static_cast<scalar>(1. / 6.)};

    /*---------------------------------------------------------------------------
     * k_{n} is always in the form of [ A(T) X B ] where A is a function of r'
     * and B is magnetic field and X symbol is for cross product. Hence dk{n}dT
     * can be represented as follows:
     *
     *  dk{n}dT =  d/dT [ A(T) X B ]  = dA(T)/dr (X) B
     *
     *  dk{n}dT = d/dT qop_n * ( T_n X B ) = qop_n ( d(T_n)/dT X B )
     *
     *  dk{n}dT = d/dT qop_n * ( T_n X B ) = qop_n ( d(T_n)/dT X B )
     *
     *  where,
     *  (X) symbol is column-wise cross product between matrix and vector,
     *  k1 = qop * T X B_first,
     *  k2 = qop * ( T + h / 2 * k1 ) X B_middle,
     *  k3 = qop * ( T + h / 2 * k2 ) X B_middle,
     *  k4 = qop * ( T + h * k3 ) X B_last.
    ---------------------------------------------------------------------------*/
    auto dk1dT = matrix_operator().template identity<3, 3>();
    auto dk2dT = matrix_operator().template identity<3, 3>();
    auto dk3dT = matrix_operator().template identity<3, 3>();
    auto dk4dT = matrix_operator().template identity<3, 3>();

    // Top-left 3x3 submatrix of Eq 3.12 of [JINST 4 P04016]
    dk1dT = sd.qop[0u] * mat_helper().column_wise_cross(dk1dT, sd.b_first);
    dk2dT = dk2dT + half_h * dk1dT;
    dk2dT = sd.qop[1u] * mat_helper().column_wise_cross(dk2dT, sd.b_middle);
    dk3dT = dk3dT + half_h * dk2dT;
    dk3dT = sd.qop[2u] * mat_helper().column_wise_cross(dk3dT, sd.b_middle);
    dk4dT = dk4dT + h * dk3dT;
    dk4dT = sd.qop[3u] * mat_helper().column_wise_cross(dk4dT, sd.b_last);

    auto dFdT = matrix_operator().template identity<3, 3>();
    auto dGdT = matrix_operator().template identity<3, 3>();
    dFdT = dFdT + h_6 * (dk1dT + dk2dT + dk3dT);
    dFdT = h * dFdT;
    dGdT = dGdT + h_6 * (dk1dT + 2.f * (dk2dT + dk3dT) + dk4dT);

    // Calculate dk{n}dqop where L is qop
    // d(k_n)d(qop_1) = d(qop_n * t_n X B_n)/d(qop_1)
    // = d(qop_n)/d(qop_1) * [ t_n (X) B_n ] + qop_n * [ d(t_n)/d(qop_1) X B_n ]
    //
    // Note that [ qop_n = qop_1 + h * d(qop_{n-1})/ds) ] as indicated by
    // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
    //
    // d(qop_n)/d(qop_1) = 1 + h * d(d(qop_{n-1})/ds)/d(qop_1))
    //
    // Let us assume that
    // d(d(qop_{n-1})/ds)/d(qop_1)) ~ d(d(qop_{n-1})/ds)/d(qop_{n-1}))
    vector3 dk1dqop = vector::cross(sd.t[0u], sd.b_first);
    vector3 dk2dqop = (1 + half_h * this->d2qopdsdqop(sd.qop[0u])) *
                          vector::cross(sd.t[1u], sd.b_middle) +
                      sd.qop[1u] * half_h * vector::cross(dk1dqop, sd.b_middle);
    vector3 dk3dqop = (1 + half_h * this->d2qopdsdqop(sd.qop[1u])) *
                          vector::cross(sd.t[2u], sd.b_middle) +
                      sd.qop[2u] * half_h * vector::cross(dk2dqop, sd.b_middle);
    vector3 dk4dqop = (1 + h * this->d2qopdsdqop(sd.qop[2u])) *
                          vector::cross(sd.t[3u], sd.b_last) +
                      sd.qop[3u] * h * vector::cross(dk3dqop, sd.b_last);

    // dFdqop and dGdqop are top-right 3x1 submatrix of equation (17)
    vector3 dFdqop = h * h_6 * (dk1dqop + dk2dqop + dk3dqop);
    vector3 dGdqop = h_6 * (dk1dqop + 2.f * (dk2dqop + dk3dqop) + dk4dqop);

    // Set transport matrix (D) and update Jacobian transport
    //( JacTransport = D * JacTransport )
    auto D = matrix_operator().template identity<e_free_size, e_free_size>();
    matrix_operator().set_block(D, dFdT, 0u, 4u);
    matrix_operator().set_block(D, dFdqop, 0u, 7u);
    matrix_operator().set_block(D, dGdT, 4u, 4u);
    matrix_operator().set_block(D, dGdqop, 4u, 7u);

    /// (4,4) element (right-bottom) of Eq. 3.12 [JINST 4 P04016]
    if (cfg.use_eloss_gradient) {
        if (this->_mat == vacuum<scalar_t>()) {
            getter::element(D, e_free_qoverp, e_free_qoverp) = 1.f;
        } else {
            // As the reference of [JINST 4 P04016] said that "The energy loss
            // and its gradient varies little within each recursion step, hence
            // the values calculated in the first stage are recycled by the
            // following stages", we obtain the d(qop)/d(qop) only from the
            // gradient at the first stage of RKN.
            //
            // But it would be better to be more precise in the future.
            getter::element(D, e_free_qoverp, e_free_qoverp) =
                1.f + d2qopdsdqop(sd.qop[0u]) * h;
        }
    }

    // Equation 3.13 of [JINST 4 P04016]
    if (cfg.use_field_gradient) {
        auto dFdR = matrix_operator().template identity<3, 3>();
        auto dGdR = matrix_operator().template identity<3, 3>();

        const vector3 pos1 = track.pos();
        const vector3 pos2 =
            pos1 + half_h * sd.t[0u] + h2 * 0.125f * sd.dtds[0u];
        const vector3 pos4 = pos1 + h * sd.t[0u] + h2 * 0.5f * sd.dtds[2u];

        const matrix_type<3, 3> field_gradient1 = evaluate_field_gradient(pos1);
        const matrix_type<3, 3> field_gradient2 = evaluate_field_gradient(pos2);
        const matrix_type<3, 3> field_gradient3 = field_gradient2;
        const matrix_type<3, 3> field_gradient4 = evaluate_field_gradient(pos4);

        // dk{n}dR = d(qop_n * t_n X B_n)/dR
        //         = qop_n * [ d(t_n)/dR (X) B_n - d(B_n)/dR (X) t_n ]
        //
        // t_n = t_1 + h * dtds_{n-1}
        // ===> d(t_n)/dR = h * dk{n-1}dR  ( Note that k == dtds )
        matrix_type<3, 3> dk1dR =
            -1.f * sd.qop[0u] *
            mat_helper().column_wise_cross(field_gradient1, sd.t[0u]);
        matrix_type<3, 3> dk2dR = sd.qop[1u] * half_h * dk1dR;
        dk2dR = mat_helper().column_wise_cross(dk2dR, sd.b_middle) -
                sd.qop[1u] *
                    mat_helper().column_wise_cross(field_gradient2, sd.t[1u]);
        matrix_type<3, 3> dk3dR = sd.qop[2u] * half_h * dk2dR;
        dk3dR = mat_helper().column_wise_cross(dk3dR, sd.b_middle) -
                sd.qop[2u] *
                    mat_helper().column_wise_cross(field_gradient3, sd.t[2u]);
        matrix_type<3, 3> dk4dR = sd.qop[3u] * h * dk3dR;
        dk4dR = mat_helper().column_wise_cross(dk4dR, sd.b_last) -
                sd.qop[3u] *
                    mat_helper().column_wise_cross(field_gradient4, sd.t[3u]);

        dFdR = dFdR + h * h_6 * (dk1dR + dk2dR + dk3dR);
        dGdR = h_6 * (dk1dR + 2.f * (dk2dR + dk3dR) + dk4dR);

        matrix_operator().set_block(D, dFdR, 0u, 0u);
        matrix_operator().set_block(D, dGdR, 4u, 0u);
    }

    this->_jac_transport = D * this->_jac_transport;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, inspector_t,
    array_t>::state::evaluate_dqopds(const std::size_t i,
                                     const typename transform3_t::scalar_type h,
                                     const scalar dqopds_prev,
                                     const detray::stepping::config& cfg) ->
    typename transform3_t::scalar_type {

    using scalar_t = typename transform3_t::scalar_type;

    const auto& track = this->_track;
    const scalar_t qop = track.qop();
    auto& sd = this->_step_data;

    if (this->_mat == detray::vacuum<scalar_t>()) {
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

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, inspector_t,
    array_t>::state::evaluate_dtds(const vector3& b_field, const std::size_t i,
                                   const typename transform3_t::scalar_type h,
                                   const vector3& dtds_prev,
                                   const typename transform3_t::scalar_type qop)
    -> vector3 {
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

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, inspector_t,
    array_t>::state::evaluate_field_gradient(const vector3& pos)
    -> matrix_type<3, 3> {

    matrix_type<3, 3> dBdr = matrix_operator().template zero<3, 3>();

    constexpr typename transform3_t::scalar_type delta =
        1e-1f * unit<scalar>::mm;

    for (unsigned int i = 0; i < 3; i++) {

        vector3 dpos1 = pos;
        dpos1[i] += delta;
        const auto bvec1_tmp =
            this->_magnetic_field.at(dpos1[0], dpos1[1], dpos1[2]);
        vector3 bvec1;
        bvec1[0u] = bvec1_tmp[0u];
        bvec1[1u] = bvec1_tmp[1u];
        bvec1[2u] = bvec1_tmp[2u];

        vector3 dpos2 = pos;
        dpos2[i] -= delta;
        const auto bvec2_tmp =
            this->_magnetic_field.at(dpos2[0], dpos2[1], dpos2[2]);
        vector3 bvec2;
        bvec2[0u] = bvec2_tmp[0u];
        bvec2[1u] = bvec2_tmp[1u];
        bvec2[2u] = bvec2_tmp[2u];

        const vector3 gradient = (bvec1 - bvec2) * (1.f / (2.f * delta));

        getter::element(dBdr, 0u, i) = gradient[0u];
        getter::element(dBdr, 1u, i) = gradient[1u];
        getter::element(dBdr, 2u, i) = gradient[2u];
    }

    return dBdr;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                        inspector_t, array_t>::state::dqopds() const ->
    typename transform3_t::scalar_type {
    return this->_step_data.dqopds[3u];
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                        inspector_t, array_t>::state::dqopds(const scalar_type
                                                                 qop) const ->
    typename transform3_t::scalar_type {

    using scalar_t = typename transform3_t::scalar_type;

    const auto& mat = this->_mat;

    // d(qop)ds is zero for empty space
    if (mat == detray::vacuum<scalar_t>()) {
        return 0.f;
    }

    const auto pdg = this->_pdg;
    const scalar_t q = this->_track.charge();
    const scalar_t p = q / qop;
    const scalar_t mass = this->_mass;
    const scalar_t E = math::sqrt(p * p + mass * mass);

    // Compute stopping power
    const scalar_t stopping_power =
        interaction<scalar_t>().compute_stopping_power(mat, pdg, mass, qop, q);

    // Assert that a momentum is a positive value
    assert(p >= 0.f);

    // d(qop)ds, which is equal to (qop) * E * (-dE/ds) / p^2
    return qop * E * stopping_power / (p * p);
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                        inspector_t, array_t>::state::dgdqop(const scalar_type
                                                                 qop) const ->
    typename transform3_t::scalar_type {

    using scalar_t = typename transform3_t::scalar_type;

    if (this->_mat == vacuum<scalar_t>()) {
        return 0.f;
    }

    auto& track = this->_track;
    const scalar_t q = track.charge();

    // dg/d(qop) = -1 * derivation of stopping power
    const scalar_t dgdqop =
        -1.f * interaction<scalar_t>().derive_stopping_power(
                   this->_mat, this->_pdg, this->_mass, qop, q);

    return dgdqop;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, inspector_t,
    array_t>::state::d2qopdsdqop(const scalar_type qop) const ->
    typename transform3_t::scalar_type {

    using scalar_t = typename transform3_t::scalar_type;

    if (this->_mat == vacuum<scalar_t>()) {
        return 0.f;
    }

    auto& track = this->_track;
    const scalar_t q = track.charge();
    const scalar_t p = q / qop;
    const scalar_t p2 = p * p;

    const auto& mass = this->_mass;
    const scalar_t E2 = p2 + mass * mass;

    // Interaction object
    interaction<scalar_t> I;

    // g = dE/ds = -1 * (-dE/ds) = -1 * stopping power
    const scalar_type g =
        -1.f * I.compute_stopping_power(this->_mat, this->_pdg, mass, qop, q);
    // dg/d(qop) = -1 * derivation of stopping power
    const scalar_t dgdqop = this->dgdqop(qop);

    // d(qop)/ds = - qop^3 * E * g / q^2
    const scalar_t dqopds = this->dqopds(qop);

    // Check Eq 3.12 of
    // (https://iopscience.iop.org/article/10.1088/1748-0221/4/04/P04016/meta)
    return dqopds * (1.f / qop * (3.f - p2 / E2) + 1.f / g * dgdqop);
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
template <typename propagation_state_t>
bool detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                        inspector_t,
                        array_t>::step(propagation_state_t& propagation,
                                       const detray::stepping::config& cfg) {

    using scalar_t = typename transform3_t::scalar_type;

    // Get stepper and navigator states
    state& stepping = propagation._stepping;
    auto& magnetic_field = stepping._magnetic_field;
    auto& navigation = propagation._navigation;

    auto vol = detector_volume{*navigation.detector(), navigation.volume()};
    stepping._mat = vol.material();

    auto& sd = stepping._step_data;

    scalar_t error_estimate{0.f};

    // First Runge-Kutta point
    const vector3 pos = stepping().pos();
    const auto bvec = magnetic_field.at(pos[0], pos[1], pos[2]);
    sd.b_first[0] = bvec[0];
    sd.b_first[1] = bvec[1];
    sd.b_first[2] = bvec[2];

    // qop should be recalcuated at every point
    // Reference: Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
    sd.dqopds[0u] = stepping.evaluate_dqopds(0u, 0.f, 0.f, cfg);
    sd.dtds[0u] = stepping.evaluate_dtds(sd.b_first, 0u, 0.f,
                                         vector3{0.f, 0.f, 0.f}, sd.qop[0u]);

    const auto try_rk4 = [&](const scalar& h) -> bool {
        // State the square and half of the step size
        const scalar_t h2{h * h};
        const scalar_t half_h{h * 0.5f};

        // Second Runge-Kutta point
        // qop should be recalcuated at every point
        // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
        const vector3 pos1 =
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
        const vector3 pos2 = pos + h * sd.t[0u] + h2 * 0.5f * sd.dtds[2u];
        const auto bvec2 = magnetic_field.at(pos2[0], pos2[1], pos2[2]);
        sd.b_last[0] = bvec2[0];
        sd.b_last[1] = bvec2[1];
        sd.b_last[2] = bvec2[2];

        sd.dqopds[3u] = stepping.evaluate_dqopds(3u, h, sd.dqopds[2u], cfg);
        sd.dtds[3u] =
            stepping.evaluate_dtds(sd.b_last, 3u, h, sd.dtds[2u], sd.qop[3u]);

        // Compute and check the local integration error estimate
        // @Todo
        const vector3 err_vec =
            h2 * (sd.dtds[0u] - sd.dtds[1u] - sd.dtds[2u] + sd.dtds[3u]);
        error_estimate =
            math::max(getter::norm(err_vec), static_cast<scalar_t>(1e-20));

        return (error_estimate <= cfg.rk_error_tol);
    };

    // Initial step size estimate
    stepping.set_step_size(navigation());

    scalar_t step_size_scaling{1.f};
    std::size_t n_step_trials{0u};

    // Adjust initial step size to integration error
    while (!try_rk4(stepping._step_size)) {

        step_size_scaling = math::min(
            math::max(0.25f * unit<scalar_t>::mm,
                      math::sqrt(math::sqrt(
                          (cfg.rk_error_tol / math::abs(error_estimate))))),
            static_cast<scalar_t>(4));

        // Only step size reduction is allowed so that we don't overstep
        assert(step_size_scaling <= 1.f);

        stepping._step_size *= step_size_scaling;

        // If step size becomes too small the particle remains at the
        // initial place
        if (math::abs(stepping._step_size) < math::abs(cfg.min_stepsize)) {
            // Not moving due to too low momentum needs an aborter
            return navigation.abort();
        }

        // If the parameter is off track too much or given step_size is not
        // appropriate
        if (n_step_trials > cfg.max_rk_updates) {
            // Too many trials, have to abort
            return navigation.abort();
        }
        n_step_trials++;

        // Run inspection while the stepsize is getting adjusted
        stepping.run_inspector(cfg, "Adjust stepsize: ", n_step_trials,
                               step_size_scaling);
    }

    // Update navigation direction
    const step::direction step_dir = stepping._step_size >= 0.f
                                         ? step::direction::e_forward
                                         : step::direction::e_backward;
    stepping.set_direction(step_dir);

    // Check constraints
    if (math::abs(stepping.step_size()) >
        math::abs(
            stepping.constraints().template size<>(stepping.direction()))) {

        // Run inspection before step size is cut
        stepping.run_inspector(cfg, "Before constraint: ");

        stepping.set_step_size(
            stepping.constraints().template size<>(stepping.direction()));
    }

    // Advance track state
    stepping.advance_track();

    // Advance jacobian transport
    stepping.advance_jacobian(cfg);

    // Call navigation update policy
    policy_t{}(stepping.policy_state(), propagation);

    // Run final inspection
    stepping.run_inspector(cfg, "Step complete: ");

    return true;
}
