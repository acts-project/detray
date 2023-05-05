/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algebra include(s).
#include "detray/plugins/algebra/array_definitions.hpp"

// Project include(s).
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/simulation/simulator.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

using namespace detray;
using transform3 = __plugin::transform3<detray::scalar>;

int main() {

    // Memory resource
    vecmem::host_memory_resource host_mr;

    // Create B field
    const vector3 B{0.f, 0.f, 2.f * unit<scalar>::T};

    // Create geometry
    using b_field_t = decltype(create_toy_geometry(host_mr))::bfield_type;
    const auto detector = create_toy_geometry(
        host_mr,
        b_field_t(b_field_t::backend_t::configuration_t{B[0], B[1], B[2]}));

    // Create track generator
    constexpr std::size_t theta_steps{4u};
    constexpr std::size_t phi_steps{4u};
    const vector3 ori{0.f, 0.f, 0.f};
    const scalar mom{1.f * unit<scalar>::GeV};
    auto generator = uniform_track_generator<free_track_parameters<transform3>>(
        theta_steps, phi_steps, ori, mom);

    // Create smearer
    measurement_smearer<transform3> smearer(100.f * unit<scalar>::um,
                                            100.f * unit<scalar>::um);

    std::size_t n_events = 2u;
    auto sim = simulator(n_events, detector, std::move(generator), smearer);

    sim.run();

    return 0;
}