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

    using scalar_t = typename transform3_t::scalar_type;

    const auto& sd = this->_step_data;
    const scalar_t h{this->_step_size};
    // const auto& mass = this->_mass;
    auto& track = this->_track;

    // Half step length
    const scalar_t half_h{h * 0.5f};
    const scalar_t h_6{h * static_cast<scalar>(1. / 6.)};

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

    // Top-left 3x3 submatrix of Eq 3.12 of [JINST 4 P04016]
    dk1dT = qop1 * mat_helper().column_wise_cross(dk1dT, sd.b_first);
    dk2dT = dk2dT + half_h * dk1dT;
    dk2dT = qop2 * mat_helper().column_wise_cross(dk2dT, sd.b_middle);
    dk3dT = dk3dT + half_h * dk2dT;
    dk3dT = qop3 * mat_helper().column_wise_cross(dk3dT, sd.b_middle);
    dk4dT = dk4dT + h * dk3dT;
    dk4dT = qop4 * mat_helper().column_wise_cross(dk4dT, sd.b_last);

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

    /// (4,4) element (right-bottom) of Eq. 3.12 [JINST 4 P04016]
    if (cfg.use_eloss_gradient) {
        if (this->_mat == vacuum<scalar_t>()) {
            getter::element(D, e_free_qoverp, e_free_qoverp) = 1.f;
        } else {
            const scalar_t q = track.charge();
            const scalar_t p = q / qop1;
            const auto& mass = this->_mass;
            const scalar_t p2 = p * p;
            const scalar_t E2 = p2 + mass * mass;
            const scalar_t E = math_ns::sqrt(E2);

            // Interaction object
            interaction<scalar> I;

            // g: dE/ds = -1 * stopping power
            const scalar_t g =
                -1.f *
                I.compute_stopping_power(this->_mat, this->_pdg, mass, qop1, q);
            // dg/d(qop) = -1 * derivation of stopping power
            const scalar_t dg_dQop =
                -1.f *
                I.derive_stopping_power(this->_mat, this->_pdg, mass, qop1, q);

            // d(Qop)/ds = - qop^3 * E * g / q^2
            const scalar_t dQop_ds =
                (-1.f * qop1 * qop1 * qop1 * E * g) / (q * q);

            const scalar_t gradient =
                dQop_ds * (1.f / qop1 * (3.f - p2 / E2) + 1.f / g * dg_dQop);

            // As the reference of [JINST 4 P04016] said that "The energy loss
            // and its gradient varies little within each recursion step, hence
            // the values calculated in the first stage are recycled by the
            // following stages", we obtain the d(qop)/d(qop) only from the
            // gradient at the first stage of RKN.
            //
            // But it would be better to be more precise in the future.
            getter::element(D, 7u, 7u) = 1.f + gradient * h;
        }
    }

    // Equation 3.13 of [JINST 4 P04016]
    if (cfg.use_field_gradient) {
        auto dFdR = matrix_operator().template identity<3, 3>();
        auto dGdR = matrix_operator().template identity<3, 3>();

        const vector3 pos1 = track.pos();
        const vector3 pos2 = pos1 + half_h * dir1;
        const vector3 pos3 = pos1 + half_h * dir2;
        const vector3 pos4 = pos1 + h * dir3;

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
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, inspector_t,
    array_t>::state::evaluate_qop(const typename transform3_t::scalar_type h,
                                  const detray::stepping::config& cfg) const ->
    typename transform3_t::scalar_type {

    using scalar_t = typename transform3_t::scalar_type;

    const auto& track = this->_track;
    const scalar_t qop = track.qop();

    if (cfg.use_mean_loss) {

        const auto& mat = this->_mat;
        if (mat == detray::vacuum<scalar_t>()) {
            return qop;
        }

        const scalar_t mass = this->_mass;
        const auto pdg = this->_pdg;
        const scalar_t p = track.p();
        const scalar_t q = track.charge();
        const auto direction = this->_direction;

        const scalar_t eloss =
            interaction<scalar_t>().compute_energy_loss_bethe(h, mat, pdg, mass,
                                                              qop, q);

        // @TODO: Recycle the codes in pointwise_material_interactor.hpp
        // Get new Energy
        const scalar_t nextE{
            std::sqrt(mass * mass + p * p) -
            std::copysign(eloss, static_cast<scalar_t>(direction))};

        // Put particle at rest if energy loss is too large
        const scalar_t nextP{
            (mass < nextE) ? std::sqrt(nextE * nextE - mass * mass) : 0.f};

        constexpr scalar_t inv{detail::invalid_value<scalar_t>()};
        return (nextP == 0.f) ? inv : (q != 0.f) ? q / nextP : 1.f / nextP;
    }

    return qop;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, inspector_t,
    array_t>::state::evaluate_k(const vector3& b_field, const int i,
                                const typename transform3_t::scalar_type h,
                                const vector3& k_prev,
                                const typename transform3_t::scalar_type k_qop)
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
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, inspector_t,
    array_t>::state::evaluate_field_gradient(const vector3& pos)
    -> matrix_type<3, 3> {

    matrix_type<3, 3> dBdr = matrix_operator().template zero<3, 3>();

    constexpr typename transform3_t::scalar_type h = 1e-4f;

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
          typename constraint_t, typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
auto detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                        inspector_t, array_t>::state::dqopds() const ->
    typename transform3_t::scalar_type {

    using scalar_t = typename transform3_t::scalar_type;

    const auto& mat = this->_mat;
    const auto pdg = this->_pdg;

    // d(qop)ds is zero for empty space
    if (mat == detray::vacuum<scalar_t>()) {
        return 0.f;
    }

    const scalar_t qop = this->_step_data.k_qop4;
    const scalar_t q = this->_track.charge();
    const scalar_t p = q / qop;
    const scalar_t mass = this->_mass;
    const scalar_t E = std::sqrt(p * p + mass * mass);

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
    const vector3 dir = stepping().dir();
    const auto bvec = magnetic_field.at(pos[0], pos[1], pos[2]);
    sd.b_first[0] = bvec[0];
    sd.b_first[1] = bvec[1];
    sd.b_first[2] = bvec[2];

    sd.k1 = stepping.evaluate_k(sd.b_first, 0, 0.f, vector3{0.f, 0.f, 0.f},
                                stepping().qop());

    const auto try_rk4 = [&](const scalar& h) -> bool {
        // State the square and half of the step size
        const scalar_t h2{h * h};
        const scalar_t half_h{h * 0.5f};

        // Second Runge-Kutta point
        const vector3 pos1 = pos + half_h * dir + h2 * 0.125f * sd.k1;
        const auto bvec1 = magnetic_field.at(pos1[0], pos1[1], pos1[2]);
        sd.b_middle[0] = bvec1[0];
        sd.b_middle[1] = bvec1[1];
        sd.b_middle[2] = bvec1[2];

        sd.k_qop2 = stepping.evaluate_qop(half_h, cfg);
        sd.k2 = stepping.evaluate_k(sd.b_middle, 1, half_h, sd.k1, sd.k_qop2);

        // Third Runge-Kutta point
        sd.k3 = stepping.evaluate_k(sd.b_middle, 2, half_h, sd.k2, sd.k_qop2);

        // Last Runge-Kutta point
        const vector3 pos2 = pos + h * dir + h2 * 0.5f * sd.k3;
        const auto bvec2 = magnetic_field.at(pos2[0], pos2[1], pos2[2]);
        sd.b_last[0] = bvec2[0];
        sd.b_last[1] = bvec2[1];
        sd.b_last[2] = bvec2[2];
        sd.k_qop4 = stepping.evaluate_qop(h, cfg);
        sd.k4 = stepping.evaluate_k(sd.b_last, 3, h, sd.k3, sd.k_qop4);

        // Compute and check the local integration error estimate
        // @Todo
        const vector3 err_vec = h2 * (sd.k1 - sd.k2 - sd.k3 + sd.k4);
        error_estimate =
            std::max(getter::norm(err_vec), static_cast<scalar_t>(1e-20));

        return (error_estimate <= cfg.rk_error_tol);
    };

    // Initial step size estimate
    stepping.set_step_size(navigation());

    scalar_t step_size_scaling{1.f};
    std::size_t n_step_trials{0u};

    // Adjust initial step size to integration error
    while (!try_rk4(stepping._step_size)) {

        step_size_scaling = std::min(
            std::max(0.25f * unit<scalar_t>::mm,
                     std::sqrt(std::sqrt(
                         (cfg.rk_error_tol / std::abs(2.f * error_estimate))))),
            static_cast<scalar_t>(4));

        // Only step size reduction is allowed so that we don't overstep
        assert(step_size_scaling <= 1.f);

        stepping._step_size *= step_size_scaling;

        // If step size becomes too small the particle remains at the
        // initial place
        if (std::abs(stepping._step_size) < std::abs(cfg.min_stepsize)) {
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
    if (std::abs(stepping.step_size()) >
        std::abs(
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
