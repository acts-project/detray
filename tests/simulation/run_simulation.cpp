/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algebra include(s).
#include "detray/plugins/algebra/array_definitions.hpp"

// Project include(s).
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/simulation/simulator.hpp"
#include "detray/simulation/track_generators.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

using namespace detray;
using transform3 = __plugin::transform3<detray::scalar>;

int main() {

    // Memory resource
    vecmem::host_memory_resource host_mr;

    // Create B field
    const vector3 B{0, 0, 2 * unit_constants::T};

    // Create geometry
    using b_field_t = decltype(create_toy_geometry(host_mr))::bfield_type;
    const auto detector = create_toy_geometry(
        host_mr,
        b_field_t(b_field_t::backend_t::configuration_t{B[0], B[1], B[2]}));

    // Create track generator
    constexpr unsigned int theta_steps{4};
    constexpr unsigned int phi_steps{4};
    const vector3 ori{0, 0, 0};
    const scalar mom = 1 * unit_constants::GeV;
    auto generator = uniform_track_generator<free_track_parameters<transform3>>(
        theta_steps, phi_steps, ori, mom);

    // Create smearer
    measurement_smearer<scalar> smearer(100 * unit_constants::um,
                                        100 * unit_constants::um);

    std::size_t n_events = 2;
    auto sim = simulator(n_events, detector, generator, smearer);

    sim.run();

    return 0;
}