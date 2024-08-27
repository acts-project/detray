/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/pdg_particle.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/constrained_step.hpp"
#include "detray/propagator/stepping_config.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/curvilinear_frame.hpp"

namespace detray {

namespace stepping {

/// A void inpector that does nothing.
///
/// Inspectors can be plugged in to understand the current stepper state.
struct void_inspector {
    template <typename state_t>
    DETRAY_HOST_DEVICE constexpr void operator()(const state_t & /*ignored*/,
                                                 const char * /*ignored*/) {
        /*Do nothing*/
    }
};

}  // namespace stepping

/// Base stepper implementation
template <typename algebra_t, typename constraint_t, typename policy_t,
          typename inspector_t = stepping::void_inspector>
class base_stepper {

    public:
    using inspector_type = inspector_t;
    using policy_type = policy_t;

    using scalar_type = dscalar<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;

    using matrix_operator = dmatrix_operator<algebra_t>;

    using free_track_parameters_type = free_track_parameters<algebra_t>;
    using bound_track_parameters_type = bound_track_parameters<algebra_t>;

    using free_matrix_type = free_matrix<algebra_t>;
    using bound_matrix_type = bound_matrix<algebra_t>;
    using bound_to_free_matrix_type = bound_to_free_matrix<algebra_t>;

    /// @brief State struct holding the track
    ///
    /// @note It has to cast into a const track via the call operation.
    struct state {

        /// Sets track parameters.
        DETRAY_HOST_DEVICE
        explicit state(const free_track_parameters_type &free_params)
            : _track(free_params) {

            curvilinear_frame<algebra_t> cf(free_params);

            // Set bound track parameters
            _bound_params.set_parameter_vector(cf.m_bound_vec);

            // A dummy covariance - should not be used
            _bound_params.set_covariance(
                matrix_operator()
                    .template identity<e_bound_size, e_bound_size>());

            // A dummy barcode - should not be used
            _bound_params.set_surface_link(geometry::barcode{});

            // Set the bound to free jacobian
            _jac_to_global = cf.bound_to_free_jacobian();

            // Reset the path length
            _s = 0.f;

            // Reset jacobian transport to identity matrix
            matrix_operator().set_identity(_jac_transport);
        }

        /// Sets track parameters from bound track parameter.
        template <typename detector_t>
        DETRAY_HOST_DEVICE state(
            const bound_track_parameters_type &bound_params,
            const detector_t &det)
            : _bound_params(bound_params) {

            // Surface
            const auto sf = tracking_surface{det, bound_params.surface_link()};

            const typename detector_t::geometry_context ctx{};
            sf.template visit_mask<
                typename parameter_resetter<algebra_t>::kernel>(
                sf.transform(ctx), *this);
        }

        /// free track parameter
        free_track_parameters_type _track;

        /// Full jacobian
        bound_matrix_type _full_jacobian =
            matrix_operator().template identity<e_bound_size, e_bound_size>();

        /// jacobian transport matrix
        free_matrix_type _jac_transport =
            matrix_operator().template identity<e_free_size, e_free_size>();

        /// bound-to-free jacobian from departure surface
        bound_to_free_matrix_type _jac_to_global =
            matrix_operator().template zero<e_free_size, e_bound_size>();

        /// bound covariance
        bound_track_parameters_type _bound_params;

        /// @returns track parameters - const access
        DETRAY_HOST_DEVICE
        free_track_parameters_type &operator()() { return _track; }

        /// @returns track parameters.
        DETRAY_HOST_DEVICE
        const free_track_parameters_type &operator()() const { return _track; }

        step::direction _direction{step::direction::e_forward};

        // Stepping constraints
        constraint_t _constraint = {};

        // Navigation policy state
        typename policy_t::state _policy_state = {};

        /// The inspector type of the stepping
        inspector_type _inspector;

        /// Track path length
        scalar_type _path_length{0.f};

        /// Absolute path length
        scalar_type _abs_path_length{0.f};

        /// Track path length from the last surface. It will be reset to 0 when
        /// the track reaches a new surface
        scalar_type _s{0.f};

        /// Current step size
        scalar_type _step_size{0.f};

        /// Previous step size (DEBUG purpose only)
        scalar_type _prev_step_size{0.f};

        /// The default particle hypothesis is muon
        pdg_particle<scalar_type> _ptc = muon<scalar_type>();

        /// is step size just initialized
        bool _initialized = true;

        /// Set new step constraint
        template <step::constraint type = step::constraint::e_actor>
        DETRAY_HOST_DEVICE inline void set_constraint(scalar_type step_size) {
            _constraint.template set<type>(step_size);
        }

        /// Set new navigation direction
        DETRAY_HOST_DEVICE inline void set_direction(step::direction dir) {
            _direction = dir;
        }

        /// @returns access to this states step constraints
        DETRAY_HOST_DEVICE
        inline const constraint_t &constraints() const { return _constraint; }

        /// @returns access to this states step constraints
        DETRAY_HOST_DEVICE
        inline typename policy_t::state &policy_state() {
            return _policy_state;
        }

        /// @returns the navigation direction
        DETRAY_HOST_DEVICE
        inline step::direction direction() const { return _direction; }

        /// Remove [all] constraints
        template <step::constraint type = step::constraint::e_actor>
        DETRAY_HOST_DEVICE inline void release_step() {
            _constraint.template release<type>();
        }

        /// Set next step size
        DETRAY_HOST_DEVICE
        inline void set_step_size(const scalar_type step) { _step_size = step; }

        /// @returns the current step size of this state.
        DETRAY_HOST_DEVICE
        inline scalar_type step_size() const { return _step_size; }

        /// @returns this states remaining path length.
        DETRAY_HOST_DEVICE
        inline scalar_type path_length() const { return _path_length; }

        /// @returns the stepping inspector
        DETRAY_HOST
        inline constexpr auto &inspector() { return _inspector; }

        /// Call the stepping inspector
        DETRAY_HOST_DEVICE
        inline void run_inspector([[maybe_unused]] const stepping::config &cfg,
                                  [[maybe_unused]] const char *message) {
            if constexpr (!std::is_same_v<inspector_t,
                                          stepping::void_inspector>) {
                _inspector(*this, cfg, message);
            }
        }
    };
};

}  // namespace detray
