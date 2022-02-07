/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/core/track.hpp"
#include "detray/tools/line_stepper.hpp"
#include "detray/tools/navigator.hpp"
#include "detray/utils/indexing.hpp"

#include "detray/core/mask_store.hpp"
#include "detray/core/transform_store.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/object_registry.hpp"
#include "detray/geometry/surface_base.hpp"
#include "detray/geometry/volume.hpp"
#include "detray/masks/masks.hpp"

#include <vecmem/memory/host_memory_resource.hpp>

/// @note __plugin has to be defined with a preprocessor command
namespace detray {

using vector3 = __plugin::vector3<detray::scalar>;

/** A navigation inspector that relays information about the encountered
 *  objects the way we need them to compare with the ray
 */
template <int navigation_status = 0,
          template <typename...> class vector_t = dvector>
struct object_tracer {

    // record all object id the navigator encounters
    vector_t<intersection> object_trace = {};

    template <typename state_type>
    auto operator()(state_type &state, const char * /*message*/) {
        // Record the candidate of an encountered object
        if (state.status() == navigation_status) {
            object_trace.push_back(std::move(*(state.current())));
        }
    }

    auto operator[](std::size_t i) { return object_trace[i]; }
};

/** A navigation inspector that prints information about the current navigation
 * state. Meant for debugging.
 */
struct print_inspector {

    // Debug output if an error in the trace is discovered
    std::stringstream debug_stream;

    template <typename state_type>
    auto operator()(state_type &state, const char *message) {
        std::string msg(message);

        debug_stream << msg << std::endl;

        debug_stream << "Volume\t\t\t\t\t\t" << state.volume() << std::endl;
        debug_stream << "surface kernel size\t\t" << state.candidates().size()
                     << std::endl;

        debug_stream << "Surface candidates: " << std::endl;
        for (const auto &sf_cand : state.candidates()) {
            debug_stream << sf_cand.to_string();
        }
        if (not state.candidates().empty()) {
            debug_stream << "=> next: ";
            if (state.is_exhausted()) {
                debug_stream << "exhausted" << std::endl;
            } else {
                debug_stream << " -> " << state.next()->index << std::endl;
            }
        }

        switch (static_cast<int>(state.status())) {
            case -3:
                debug_stream << "status\t\t\t\t\ton_target" << std::endl;
                break;
            case -2:
                debug_stream << "status\t\t\t\t\tabort" << std::endl;
                break;
            case -1:
                debug_stream << "status\t\t\t\t\tunknowm" << std::endl;
                break;
            case 0:
                debug_stream << "status\t\t\t\t\ttowards_surface" << std::endl;
                break;
            case 1:
                debug_stream << "status\t\t\t\t\ton_surface" << std::endl;
                break;
            case 2:
                debug_stream << "status\t\t\t\t\ttowards_portal" << std::endl;
                break;
            case 3:
                debug_stream << "status\t\t\t\t\ton_portal" << std::endl;
                break;
        };
        debug_stream << "current object\t\t" << state.on_object() << std::endl;
        debug_stream << "distance to next\t";
        if (std::abs(state()) < state.tolerance()) {
            debug_stream << "on obj (within tol)" << std::endl;
        } else {
            debug_stream << state() << std::endl;
        }
        switch (state.trust_level()) {
            case 0:
                debug_stream << "trust\t\t\t\t\tno_trust" << std::endl;
                break;
            case 1:
                debug_stream << "trust\t\t\t\t\tfair_trust" << std::endl;
                break;
            case 3:
                debug_stream << "trust\t\t\t\t\thigh_trust" << std::endl;
                break;
            case 4:
                debug_stream << "trust\t\t\t\t\tfull_trust" << std::endl;
                break;
        };
        debug_stream << std::endl;
    }

    std::string to_string() { return debug_stream.str(); }
};

/** A navigation inspector that aggregates a number of different inspectors.*/
template <typename... Inspectors>
struct aggregate_inspector {

    using inspector_tuple_t = std::tuple<Inspectors...>;
    inspector_tuple_t _inspectors{};

    template <unsigned int current_id = 0, typename state_type>
    auto operator()(state_type &state, const char *message) {
        // Call inspector
        std::get<current_id>(_inspectors)(state, message);

        // Next mask type
        if constexpr (current_id <
                      std::tuple_size<inspector_tuple_t>::value - 1) {
            return operator()<current_id + 1>(state, message);
        }
    }

    template <typename inspector_t>
    decltype(auto) get() {
        return std::get<inspector_t>(_inspectors);
    }
};

template<template <typename, unsigned int> class array_t = darray,
         template <typename...> class tuple_t = dtuple,
         template <typename...> class vector_type = dvector>
struct surface_vector {
    template <typename T>
    using vector_t = vector_type<T>;

    using transform_store = static_transform_store<vector_type>;
    using context_type = typename transform_store::context;

    using rectangle =
        rectangle2<planar_intersector, __plugin::cartesian2<detray::scalar>,
                   array_t<dindex, 2>, 0>;

    using surface_type = surface_base<dindex, array_t<dindex, 2>, dindex,
                                      dindex, array_t<dindex, 2>>;
    using volume_type = volume<object_registry<surface_type>, dindex_range, array_t>;

    surface_vector(vecmem::memory_resource &resource, dvector<scalar> distances, vector3 direction, context_type ctx = {}) 
        : _surfaces(&resource),
          _transforms(resource),
          _masks(resource)
    {
        // Rotation matrix
        vector3 z = direction;
        vector3 x = vector::normalize(vector3{0, -z[2], z[1]});
        //vector3 z = vector::normalize(vector3{0., 0., 1.});
        //vector3 x = vector::normalize(vector3{1., 0., 0.});

        const array_t<scalar, 6> bounds = {};
        volume_type &vol = _volumes.emplace_back(bounds);
        vol.set_index(_volumes.size() - 1);
        _surfaces.reserve(distances.size());
        _masks.template add_mask<0>(10., 10.);
        array_t<dindex, 2> mask_link{0, _masks.template size<0>()};
        for (auto &d : distances) {
            vector3 t = d * direction;
            _transforms.emplace_back(ctx, t, z, x);
            _surfaces.emplace_back(_transforms.size(ctx) - 1, mask_link, 0, dindex_invalid);
            _surfaces.back().set_edge({0, 0});
        }
        // last surface is volume portal
        _surfaces.back().set_edge({0, dindex_invalid});
        vol.update_range({0, _surfaces.size()});
    }
    // navigator interface
    inline auto &volumes() const { return _volumes; }
    inline auto &surfaces() const { return _surfaces; }
    inline const auto &transforms(const context_type & /*ctx*/ = {}) const {
        return _transforms;
    }
    inline auto &masks() const { return _masks; }

    vector_t<volume_type> _volumes = {};
    vector_t<surface_type> _surfaces = {};
    transform_store _transforms = {};
    mask_store<tuple_t, vector_t, rectangle> _masks;
};

} // namespace detray

// This tests the construction and general methods of the navigator
TEST(ALGEBRA_PLUGIN, guided_navigator) {
    using namespace detray;
    vecmem::host_memory_resource host_mr;

    dvector<scalar> dists = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
    surface_vector<> surfaces(host_mr, dists, vector::normalize(vector3{1., 0., 0.}));
    
    using detray_inspector =
        aggregate_inspector<object_tracer<1>, print_inspector>;
    using guided_navigator = navigator<decltype(surfaces), detray_inspector>;
    using nav_context = decltype(surfaces)::context_type;
    using stepper = line_stepper<track<nav_context>>;

    guided_navigator n(surfaces);
    stepper s;

    // test track
    track<nav_context> traj;
    traj.pos = {0., 0., 0.};
    traj.dir = vector::normalize(vector3{1., 0., 0.});
    traj.ctx = nav_context{};
    traj.momentum = 100.;
    traj.overstep_tolerance = -1e-4;

    stepper::state s_state(traj);
    guided_navigator::state n_state;

    // Set initial volume (no grid yet)
    n_state.set_volume(0u);

    bool heartbeat = n.status(n_state, traj);
    // Run while there is a heartbeat
    while (heartbeat) {
        // (Re-)target
        heartbeat &= n.target(n_state, s_state());
        // Take the step
        heartbeat &= s.step(s_state, n_state());
        // And check the status
        heartbeat &= n.status(n_state, s_state());
    }

    // sequence of surfaces we expect to see
    std::vector<dindex> sf_sequence = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    auto &obj_tracer =
        n_state.inspector().template get<object_tracer<1>>();
    auto &debug_printer =
        n_state.inspector().template get<print_inspector>();
    EXPECT_EQ(obj_tracer.object_trace.size(), sf_sequence.size())<< debug_printer.to_string();

    // Every iteration steps through one barrel layer
    for (const auto &sf_id : sf_sequence) {
    }
}
