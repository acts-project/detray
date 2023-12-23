/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// System include(s)
#include <cmath>

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t,
          template <typename, std::size_t> class array_t>
void detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                        array_t>::state::advance_track() {

    const auto& sd = this->_step_data;
    const scalar h{this->_step_size};
    const scalar h_6{h * static_cast<scalar>(1. / 6.)};
    auto& track = this->_track;
    auto pos = track.pos();
    auto dir = track.dir();

    // Update the track parameters according to the equations of motion
    pos = pos + h * (dir + h_6 * (sd.k1 + sd.k2 + sd.k3));
    track.set_pos(pos);

    dir = dir + h_6 * (sd.k1 + 2.f * (sd.k2 + sd.k3) + sd.k4);
    dir = vector::normalize(dir);
    track.set_dir(dir);
    track.set_qop(sd.k_qop4);

    // Update path length
    this->_path_length += h;
    this->_s += h;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t,
          template <typename, std::size_t> class array_t>
void detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                        array_t>::state::advance_jacobian() {
    /// The calculations are based on ATL-SOFT-PUB-2009-002. The update of the
    /// Jacobian matrix is requires only the calculation of eq. 17 and 18.
    /// Since the terms of eq. 18 are currently 0, this matrix is not needed
    /// in the calculation. The matrix A from eq. 17 consists out of 3
    /// different parts. The first one is given by the upper left 3x3 matrix
    /// that are calculated by the derivatives dF/dT (called dFdT) and dG/dT
    /// (calles dGdT). The second is given by the top 3 lines of the rightmost
    /// column. This is calculated by dFdL and dGdL. The remaining non-zero term
    /// is calculated directly. The naming of the variables is explained in eq.
    /// 11 and are directly related to the initial problem in eq. 7.
    /// The evaluation is based by propagating the parameters T and lambda as
    /// given in eq. 16 and evaluating the derivations for matrix A.
    /// @note The translation for u_{n+1} in eq. 7 is in this case a
    /// 3-dimensional vector without a dependency of Lambda or lambda neither in
    /// u_n nor in u_n'. The second and fourth eq. in eq. 14 have the constant
    /// offset matrices h * Id and Id respectively. This involves that the
    /// constant offset does not exist for rectangular matrix dGdu' (due to the
    /// missing Lambda part) and only exists for dFdu' in dlambda/dlambda.

    const auto& sd = this->_step_data;
    const scalar h{this->_step_size};
    // const auto& mass = this->_mass;
    auto& track = this->_track;

    // Half step length
    const scalar half_h{h * 0.5f};
    const scalar h_6{h * static_cast<scalar>(1. / 6.)};

    // Direction
    const auto dir1 = track.dir();
    const vector3 dir2 = dir1 + half_h * sd.k1;
    const vector3 dir3 = dir1 + half_h * sd.k2;
    const vector3 dir4 = dir1 + h * sd.k3;

    // Q over P
    const auto qop1 = track.qop();
    const auto qop2 = sd.k_qop2;
    const auto qop3 = qop2;
    const auto qop4 = sd.k_qop4;

    /*---------------------------------------------------------------------------
     * k_{n} is always in the form of [ A(T) X B ] where A is a function of r'
     * and B is magnetic field and X symbol is for cross product. Hence dk{n}dT
     * can be represented as follows:
     *
     *  dk{n}dT =  d/dT [ A(T) X B ] = dA(T)/dr (X) B
     *
     *  where,
     *  (X) symbol is column-wise cross product between matrix and vector,
     *  k1 = qop * T X B_first,
     *  k2 = qop * ( T + h / 2 * k1 ) X B_middle,
     *  k3 = qop * ( T + h / 2 * k2 ) X B_middle,
     *  k4 = qop * ( T + h * k3 ) * B_last.
    ---------------------------------------------------------------------------*/
    auto dk1dT = matrix_operator().template identity<3, 3>();
    auto dk2dT = matrix_operator().template identity<3, 3>();
    auto dk3dT = matrix_operator().template identity<3, 3>();
    auto dk4dT = matrix_operator().template identity<3, 3>();

    dk1dT = qop1 * mat_helper().column_wise_cross(dk1dT, sd.b_first);
    dk2dT = dk2dT + half_h * dk1dT;
    dk2dT = qop2 * mat_helper().column_wise_cross(dk2dT, sd.b_middle);
    dk3dT = dk3dT + half_h * dk2dT;
    dk3dT = qop3 * mat_helper().column_wise_cross(dk3dT, sd.b_middle);
    dk4dT = dk4dT + h * dk3dT;
    dk4dT = qop4 * mat_helper().column_wise_cross(dk4dT, sd.b_last);

    // dFdT and dGdT are top-left 3x3 submatrix of equation (17)
    auto dFdT = matrix_operator().template identity<3, 3>();
    auto dGdT = matrix_operator().template identity<3, 3>();
    dFdT = dFdT + h_6 * (dk1dT + dk2dT + dk3dT);
    dFdT = h * dFdT;
    dGdT = dGdT + h_6 * (dk1dT + 2.f * (dk2dT + dk3dT) + dk4dT);

    // Calculate dk{n}dL where L is qop
    vector3 dk1dL = vector::cross(dir1, sd.b_first);
    vector3 dk2dL = vector::cross(dir2, sd.b_middle) +
                    qop2 * half_h * vector::cross(dk1dL, sd.b_middle);
    vector3 dk3dL = vector::cross(dir3, sd.b_middle) +
                    qop3 * half_h * vector::cross(dk2dL, sd.b_middle);
    vector3 dk4dL = vector::cross(dir4, sd.b_last) +
                    qop4 * h * vector::cross(dk3dL, sd.b_last);

    // dFdL and dGdL are top-right 3x1 submatrix of equation (17)
    vector3 dFdL = h * h_6 * (dk1dL + dk2dL + dk3dL);
    vector3 dGdL = h_6 * (dk1dL + 2.f * (dk2dL + dk3dL) + dk4dL);

    // Set transport matrix (D) and update Jacobian transport
    //( JacTransport = D * JacTransport )
    auto D = matrix_operator().template identity<e_free_size, e_free_size>();
    matrix_operator().set_block(D, dFdT, 0u, 4u);
    matrix_operator().set_block(D, dFdL, 0u, 7u);
    matrix_operator().set_block(D, dGdT, 4u, 4u);
    matrix_operator().set_block(D, dGdL, 4u, 7u);

    if (_use_field_gradient) {
        auto dFdR = matrix_operator().template identity<3, 3>();
        auto dGdR = matrix_operator().template identity<3, 3>();

        const vector3 pos1 = track.pos();
        const vector3 pos2 = pos1 + half_h * dir2;
        const vector3 pos3 = pos1 + half_h * dir3;
        const vector3 pos4 = pos1 + h * dir4;

        const matrix_type<3, 3> field_gradient1 = evaluate_field_gradient(pos1);
        const matrix_type<3, 3> field_gradient2 = evaluate_field_gradient(pos2);
        const matrix_type<3, 3> field_gradient3 = evaluate_field_gradient(pos3);
        const matrix_type<3, 3> field_gradient4 = evaluate_field_gradient(pos4);

        // dk{n}dR = d(qop_n * t_n X B_n)/dR
        //         = qop_n * [ d(t_n)/dR (X) B_n - d(B_n)/dR (X) t_n ]
        matrix_type<3, 3> dk1dR =
            -1.f * qop1 * mat_helper().column_wise_cross(field_gradient1, dir1);
        matrix_type<3, 3> dk2dR = qop2 * half_h * dk1dR;
        dk2dR = mat_helper().column_wise_cross(dk2dR, sd.b_middle) -
                qop2 * mat_helper().column_wise_cross(field_gradient2, dir2);
        matrix_type<3, 3> dk3dR = qop3 * half_h * dk2dR;
        dk3dR = mat_helper().column_wise_cross(dk3dR, sd.b_middle) -
                qop3 * mat_helper().column_wise_cross(field_gradient3, dir3);
        matrix_type<3, 3> dk4dR = qop4 * h * dk3dR;
        dk4dR = mat_helper().column_wise_cross(dk4dR, sd.b_last) -
                qop4 * mat_helper().column_wise_cross(field_gradient4, dir4);

        dFdR = dFdR + h * h_6 * (dk1dR + dk2dR + dk3dR);
        dGdR = h_6 * (dk1dR + 2.f * (dk2dR + dk3dR) + dk4dR);

        matrix_operator().set_block(D, dFdR, 0u, 0u);
        matrix_operator().set_block(D, dGdR, 4u, 0u);
    }

    /// Calculate (4,4) element of equation (17)
    /// NOTE: Let's skip this element for the moment
    /// const auto p = getter::norm(track.mom());
    /// matrix_operator().element(D, 3, 7) =
    /// h * mass * mass * qop * getter::perp(vector2{1, mass / p});

    this->_jac_transport = D * this->_jac_transport;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                        array_t>::state::evaluate_qop(const scalar h)
    -> scalar {

    const auto mat = this->_mat;
    auto& track = this->_track;
    const scalar qop = track.qop();

    if (mat == detray::vacuum<scalar>()) {
        return qop;
    }

    const scalar mass = this->_mass;
    const auto pdg = this->_pdg;
    const scalar p = track.p();
    const scalar q = track.charge();
    const auto direction = this->_direction;

    scalar new_qop = track.qop();

    if (this->_use_mean_loss) {

        const scalar eloss = interaction<scalar>().compute_energy_loss_bethe(
            h, mat, pdg, mass, qop, q);

        // @TODO: Recycle the codes in pointwise_material_interactor.hpp
        // Get new Energy
        const scalar nextE{
            std::sqrt(mass * mass + p * p) -
            std::copysign(eloss, static_cast<scalar>(direction))};

        // Put particle at rest if energy loss is too large
        const scalar nextP{
            (mass < nextE) ? std::sqrt(nextE * nextE - mass * mass) : 0.f};
        new_qop = (q != 0.f) ? q / nextP : 1.f / nextP;
    }

    return new_qop;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                        array_t>::state::evaluate_k(const vector3& b_field,
                                                    const int i, const scalar h,
                                                    const vector3& k_prev,
                                                    const scalar k_qop)
    -> vector3 {
    auto& track = this->_track;
    const auto dir = track.dir();

    vector3 k_new;

    if (i == 0) {
        k_new = k_qop * vector::cross(dir, b_field);
    } else {
        k_new = k_qop * vector::cross(dir + h * k_prev, b_field);
    }

    return k_new;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                        array_t>::state::evaluate_field_gradient(const vector3&
                                                                     pos)
    -> matrix_type<3, 3> {

    matrix_type<3, 3> dBdr = matrix_operator().template zero<3, 3>();

    const scalar h = 0.001f;

    for (unsigned int i = 0; i < 3; i++) {

        vector3 dpos1 = pos;
        dpos1[i] += h;
        const auto bvec1_tmp =
            this->_magnetic_field.at(dpos1[0], dpos1[1], dpos1[2]);
        vector3 bvec1;
        bvec1[0u] = bvec1_tmp[0u];
        bvec1[1u] = bvec1_tmp[1u];
        bvec1[2u] = bvec1_tmp[2u];

        vector3 dpos2 = pos;
        dpos2[i] -= h;
        const auto bvec2_tmp =
            this->_magnetic_field.at(dpos2[0], dpos2[1], dpos2[2]);
        vector3 bvec2;
        bvec2[0u] = bvec2_tmp[0u];
        bvec2[1u] = bvec2_tmp[1u];
        bvec2[2u] = bvec2_tmp[2u];

        const vector3 gradient = (bvec1 - bvec2) * (1.f / (2.f * h));

        getter::element(dBdr, 0u, i) = gradient[0u];
        getter::element(dBdr, 1u, i) = gradient[1u];
        getter::element(dBdr, 2u, i) = gradient[2u];
    }

    return dBdr;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t,
          template <typename, std::size_t> class array_t>
template <typename propagation_state_t>
bool detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                        array_t>::step(propagation_state_t& propagation) {

    // Get stepper and navigator states
    state& stepping = propagation._stepping;
    auto& magnetic_field = stepping._magnetic_field;
    auto& navigation = propagation._navigation;

    const auto det = navigation.detector();
    stepping._mat = det->volume_by_index(navigation.volume()).material();

    auto& sd = stepping._step_data;

    scalar error_estimate{0.f};

    // First Runge-Kutta point
    const vector3 pos = stepping().pos();
    const vector3 dir = stepping().dir();
    const auto bvec = magnetic_field.at(pos[0], pos[1], pos[2]);
    sd.b_first[0] = bvec[0];
    sd.b_first[1] = bvec[1];
    sd.b_first[2] = bvec[2];

    sd.k1 = stepping.evaluate_k(sd.b_first, 0, 0.f, vector3{0.f, 0.f, 0.f},
                                stepping().qop());

    const auto try_rk4 = [&](const scalar& h) -> bool {
        // State the square and half of the step size
        const scalar h2{h * h};
        const scalar half_h{h * 0.5f};

        // Second Runge-Kutta point
        const vector3 pos1 = pos + half_h * dir + h2 * 0.125f * sd.k1;
        const auto bvec1 = magnetic_field.at(pos1[0], pos1[1], pos1[2]);
        sd.b_middle[0] = bvec1[0];
        sd.b_middle[1] = bvec1[1];
        sd.b_middle[2] = bvec1[2];

        sd.k_qop2 = stepping.evaluate_qop(half_h);
        sd.k2 = stepping.evaluate_k(sd.b_middle, 1, half_h, sd.k1, sd.k_qop2);

        // Third Runge-Kutta point
        sd.k3 = stepping.evaluate_k(sd.b_middle, 2, half_h, sd.k2, sd.k_qop2);

        // Last Runge-Kutta point
        const vector3 pos2 = pos + h * dir + h2 * 0.5f * sd.k3;
        const auto bvec2 = magnetic_field.at(pos2[0], pos2[1], pos2[2]);
        sd.b_last[0] = bvec2[0];
        sd.b_last[1] = bvec2[1];
        sd.b_last[2] = bvec2[2];
        sd.k_qop4 = stepping.evaluate_qop(h);
        sd.k4 = stepping.evaluate_k(sd.b_last, 3, h, sd.k3, sd.k_qop4);

        // Compute and check the local integration error estimate
        // @Todo
        const vector3 err_vec = h2 * (sd.k1 - sd.k2 - sd.k3 + sd.k4);
        error_estimate =
            std::max(getter::norm(err_vec), static_cast<scalar>(1e-20));

        return (error_estimate <= stepping._tolerance);
    };

    // Initial step size estimate
    stepping.set_step_size(navigation());

    scalar step_size_scaling{1.f};
    std::size_t n_step_trials{0u};

    // Adjust initial step size to integration error
    while (!try_rk4(stepping._step_size)) {

        step_size_scaling = std::min(
            std::max(0.25f * unit<scalar>::mm,
                     std::sqrt(std::sqrt((stepping._tolerance /
                                          std::abs(2.f * error_estimate))))),
            static_cast<scalar>(4));

        stepping._step_size *= step_size_scaling;

        // If step size becomes too small the particle remains at the
        // initial place
        if (std::abs(stepping._step_size) <
            std::abs(stepping._step_size_cutoff)) {
            // Not moving due to too low momentum needs an aborter
            return navigation.abort();
        }

        // If the parameter is off track too much or given step_size is not
        // appropriate
        if (n_step_trials > stepping._max_rk_step_trials) {
            // Too many trials, have to abort
            return navigation.abort();
        }
        n_step_trials++;
    }

    // Update navigation direction
    const step::direction step_dir = stepping._step_size >= 0.f
                                         ? step::direction::e_forward
                                         : step::direction::e_backward;
    stepping.set_direction(step_dir);

    // Check constraints
    if (std::abs(stepping.step_size()) >
        std::abs(
            stepping.constraints().template size<>(stepping.direction()))) {
        stepping.set_step_size(
            stepping.constraints().template size<>(stepping.direction()));
    }

    // Advance track state
    stepping.advance_track();

    // Advance jacobian transport
    stepping.advance_jacobian();

    // Call navigation update policy
    policy_t{}(stepping.policy_state(), propagation);

    return true;
}
