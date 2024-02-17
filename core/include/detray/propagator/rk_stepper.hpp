/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/interaction.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/navigation/policies.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/matrix_helper.hpp"

namespace detray {

/// Runge-Kutta-Nystrom 4th order stepper implementation
///
/// @tparam magnetic_field_t the type of magnetic field
/// @tparam track_t the type of track that is being advanced by the stepper
/// @tparam constraint_ the type of constraints on the stepper
template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t = unconstrained_step,
          typename policy_t = stepper_rk_policy,
          typename random_device_t = stepping::void_random_device,
          typename inspector_t = stepping::void_inspector,
          template <typename, std::size_t> class array_t = darray>
class rk_stepper final
    : public base_stepper<transform3_t, constraint_t, policy_t, inspector_t> {

    public:
    using base_type = base_stepper<transform3_t, constraint_t, policy_t,
                                   random_device_t, inspector_t>;

    using transform3_type = transform3_t;
    using scalar_type = typename transform3_type::scalar_type;
    using point3 = typename transform3_type::point3;
    using vector2 = typename transform3_type::point2;
    using vector3 = typename transform3_type::vector3;
    using matrix_operator = typename base_type::matrix_operator;
    using mat_helper = matrix_helper<matrix_operator>;
    using free_track_parameters_type =
        typename base_type::free_track_parameters_type;
    using bound_track_parameters_type =
        typename base_type::bound_track_parameters_type;
    using magnetic_field_type = magnetic_field_t;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    using interaction_type = interaction<scalar_type>;

    DETRAY_HOST_DEVICE
    rk_stepper() {}

    struct state : public base_type::state {

        static constexpr const stepping::id id = stepping::id::e_rk;

        DETRAY_HOST_DEVICE
        state(const free_track_parameters_type& t,
              const magnetic_field_t& mag_field)
            : base_type::state(t), _magnetic_field(mag_field) {}

        template <typename detector_t>
        DETRAY_HOST_DEVICE state(
            const bound_track_parameters_type& bound_params,
            const magnetic_field_t& mag_field, const detector_t& det)
            : base_type::state(bound_params, det), _magnetic_field(mag_field) {}

        /// stepping data required for RKN4
        struct {
            vector3 b_first{0.f, 0.f, 0.f};
            vector3 b_middle{0.f, 0.f, 0.f};
            vector3 b_last{0.f, 0.f, 0.f};
            // t = tangential direction = dr/ds
            std::array<vector3, 4u> t;
            // q/p
            std::array<scalar_type, 4u> qop;
            // dt/ds = d^2r/ds^2 = q/p ( t X B )
            std::array<vector3, 4u> dtds;
            // d(q/p)/ds
            std::array<scalar_type, 4u> dqopds;
        } _step_data;

        /// Magnetic field view
        const magnetic_field_t _magnetic_field;

        /// Material that track is passing through. Usually a volume material
        detray::material<scalar_type> _mat = detray::vacuum<scalar_type>();

        /// Update the track state by Runge-Kutta-Nystrom integration.
        DETRAY_HOST_DEVICE
        inline void advance_track(
            const stepping::config<scalar_type>& cfg = {});

        /// Update the jacobian transport from free propagation
        DETRAY_HOST_DEVICE
        inline void advance_jacobian(
            const stepping::config<scalar_type>& cfg = {});

        /// Update the covariance from multiple scattering
        DETRAY_HOST_DEVICE
        inline void add_multiple_scattering_covariance();

        /// evaulate dqopds for a given step size and material
        DETRAY_HOST_DEVICE
        inline scalar_type evaluate_dqopds(
            const std::size_t i, const typename transform3_t::scalar_type h,
            const scalar dqopds_prev);

        /// evaulate dtds for runge kutta stepping
        DETRAY_HOST_DEVICE
        inline vector3 evaluate_dtds(const vector3& b_field,
                                     const std::size_t i, const scalar_type h,
                                     const vector3& dtds_prev,
                                     const scalar_type qop);

        DETRAY_HOST_DEVICE
        inline matrix_type<3, 3> evaluate_field_gradient(const vector3& pos);

        /// Evaluate dtds, where t is the unit tangential direction
        DETRAY_HOST_DEVICE
        inline vector3 dtds() const { return this->_step_data.dtds[3u]; }

        /// Evaulate d(qop)/ds
        DETRAY_HOST_DEVICE
        inline scalar_type dqopds() const;

        DETRAY_HOST_DEVICE
        inline scalar_type dqopds(const scalar_type qop) const;

        /// Evaulate d(d(qop)/ds)dqop
        DETRAY_HOST_DEVICE
        inline scalar_type d2qopdsdqop(const scalar_type qop) const;

        /// Call the stepping inspector
        template <typename... Args>
        DETRAY_HOST_DEVICE inline void run_inspector(
            [[maybe_unused]] const stepping::config<scalar_type>& cfg,
            [[maybe_unused]] const char* message,
            [[maybe_unused]] Args&&... args) {
            if constexpr (not std::is_same_v<inspector_t,
                                             stepping::void_inspector>) {
                this->_inspector(*this, cfg, message,
                                 std::forward<Args>(args)...);
            }
        }
    };

    /// Take a step, using an adaptive Runge-Kutta algorithm.
    ///
    /// @return returning the heartbeat, indicating if the stepping is alive
    template <typename propagation_state_t>
    DETRAY_HOST_DEVICE bool step(propagation_state_t& propagation,
                                 const stepping::config<scalar_type>& cfg = {});
};

}  // namespace detray

#include "detray/propagator/rk_stepper.ipp"
