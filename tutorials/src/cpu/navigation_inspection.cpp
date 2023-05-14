/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/intersection/detail/trajectories.hpp"  // helix
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/inspectors.hpp"
#include "tests/common/tools/particle_gun.hpp"

// Example linear algebra plugin: std::array
#include "detray/tutorial/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <sstream>

/// Run the navigation through the toy detector including a detailed record
/// of the encountered surfaces (using the navigation inspectors)
int main() {
    // Toy detector
    using toy_detector_t = detray::detector<detray::toy_metadata<>>;
    using b_field_t = typename toy_detector_t::bfield_type;

    /// Type that holds the intersection information
    using intersection_t =
        detray::intersection2D<typename toy_detector_t::surface_type,
                               detray::tutorial::transform3>;

    /// Inspector that records all encountered surfaces
    using object_tracer_t = detray::navigation::object_tracer<
        intersection_t, detray::dvector,
        detray::navigation::status::e_on_module,
        detray::navigation::status::e_on_portal>;

    /// Insoector that print the navigator state from within the navigator's
    /// method calls
    using nav_print_inspector_t = detray::navigation::print_inspector;

    /// Aggregation of multiple inspectors
    using inspector_t =
        detray::aggregate_inspector<object_tracer_t, nav_print_inspector_t>;

    // Navigation with inspection
    using navigator_t =
        detray::navigator<toy_detector_t, inspector_t, intersection_t>;
    // Runge-Kutta stepper
    using stepper_t =
        detray::rk_stepper<b_field_t::view_t, detray::tutorial::transform3>;
    // Propagator with empty actor chain
    using propagator_t =
        detray::propagator<stepper_t, navigator_t, detray::actor_chain<>>;

    // Build toy detector with constant B-field in z-direction
    const detray::tutorial::vector3 B{0. * detray::unit<detray::scalar>::T,
                                      0. * detray::unit<detray::scalar>::T,
                                      2. * detray::unit<detray::scalar>::T};

    vecmem::host_memory_resource host_mr;

    const toy_detector_t det = detray::create_toy_geometry(
        host_mr,
        b_field_t(b_field_t::backend_t::configuration_t{B[0], B[1], B[2]}));

    // Build the propagator
    propagator_t prop(stepper_t{}, navigator_t{});

    // Track generation config
    // Trivial example: Single track escapes through beampipe
    constexpr std::size_t theta_steps{1};
    constexpr std::size_t phi_steps{1};
    const detray::tutorial::point3 origin{0., 0., 0.};
    constexpr detray::scalar p_mag{10. * detray::unit<detray::scalar>::GeV};

    // Iterate through uniformly distributed momentum directions
    for (auto track : detray::uniform_track_generator<
             detray::free_track_parameters<detray::tutorial::transform3>>(
             theta_steps, phi_steps, origin, p_mag)) {
        track.set_overstep_tolerance(-7. * detray::unit<detray::scalar>::um);

        // Get ground-truth helix from track
        detray::detail::helix helix(track, &B);

        // Shoot helix through the detector and record all surface intersections
        const auto intersection_trace =
            detray::particle_gun::shoot_particle(det, helix);

        // Now follow that helix with the same track and check, if we find
        // the same volumes and distances along the way
        propagator_t::state propagation(track, det.get_bfield(), det);

        // Retrieve navigation information
        auto &inspector = propagation._navigation.inspector();
        auto &obj_tracer = inspector.template get<object_tracer_t>();
        auto &debug_printer = inspector.template get<nav_print_inspector_t>();

        // Run the actual propagation
        prop.propagate(propagation);
        std::cout << debug_printer.to_string();

        // Compare helix trace to object tracer
        std::stringstream debug_stream;
        for (std::size_t intr_idx = 0; intr_idx < intersection_trace.size();
             ++intr_idx) {
            debug_stream << "-------Intersection trace\n"
                         << "helix gun: "
                         << "found in vol: "
                         << intersection_trace[intr_idx].first << ",\n\t"
                         << intersection_trace[intr_idx].second;
            debug_stream << "\nnavig.:\t" << obj_tracer[intr_idx];
        }
        std::cout << debug_stream.str() << std::endl;

        // Compare intersection records
        //
        // [...]
        //
    }
}
