/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/base_stepper.hpp"

namespace detray {

struct bound_state_updater : actor {

    using matrix_operator = standard_matrix_operator<scalar>;
    using covariance_engine = detail::covariance_engine<scalar>;
    using vector_engine = covariance_engine::vector_engine;

    struct state {};

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state &t_state,
                                       propagator_state_t &propagation) const {

        auto &navigation = propagation._navigation;

        if (navigation.is_on_module()) {

            auto &stepping = propagation._stepping;

            typename propagator_state_t::context_type ctx{};

            const auto &det = navigation.detector();
            const auto &sf_idx = navigation.current_object();
            const auto &sf = det->surface_by_index(sf_idx);
            const auto &trf_store = det->transform_store(ctx);
            const auto &trf3 =
                trf_store.contextual_transform(ctx, sf.transform());

            // Get free vector
            const auto &fvector = stepping().vector();

            // Update bound vector
            const auto bvector =
                vector_engine().free_to_bound_vector(trf3, fvector);

            stepping._bound_params.set_vector(bvector);

            // Update bound covariance
            const auto bcov = covariance_engine().bound_to_bound_covariance(
                trf3, stepping._bound_params.covariance(), fvector,
                stepping._jac_to_global, stepping._jac_transport,
                stepping._derivative);

            stepping._bound_params.set_covariance(bcov);

            // For debugging purpose
            stepping._prev_jac_transport = stepping._jac_transport;

            // Reset the jacobians in RK stepper state
            covariance_engine().reinitialize_jacobians(
                trf3, bvector, stepping._jac_to_global, stepping._jac_transport,
                stepping._derivative);
        }
    }
};

struct free_state_updater : actor {

    using matrix_operator = standard_matrix_operator<scalar>;
    using covariance_engine = detail::covariance_engine<scalar>;
    using jacobian_engine = covariance_engine::jacobian_engine;
    using vector_engine = covariance_engine::vector_engine;

    struct state {};

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state &t_state,
                                       propagator_state_t &propagation) const {

        auto &navigation = propagation._navigation;

        if (navigation.is_on_module()) {

            auto &stepping = propagation._stepping;

            typename propagator_state_t::context_type ctx{};

            const auto &det = navigation.detector();
            const auto &sf_idx = navigation.current_object();
            const auto &sf = det->surface_by_index(sf_idx);
            const auto &trf_store = det->transform_store(ctx);
            const auto &trf3 =
                trf_store.contextual_transform(ctx, sf.transform());

            // Update free vector
            const auto &bvector = stepping._bound_params.vector();
            const auto fvector =
                vector_engine().bound_to_free_vector(trf3, bvector);

            stepping().set_vector(fvector);

            // Update free cov
            const auto &bcov = stepping._bound_params.covariance();
            const auto bound_to_free_jacobian =
                jacobian_engine().bound_to_free_coordinate(trf3, fvector);

            const auto fcov = bound_to_free_jacobian * bcov *
                              bound_to_free_jacobian.transpose();

            stepping().set_covariance(fcov);
        }
    }
};

template <typename composite_actor_type>
using surface_actor = actor_chain<dtuple, bound_state_updater,
                                  composite_actor_type, free_state_updater>;

}  // namespace detray