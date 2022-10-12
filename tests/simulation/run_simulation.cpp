/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algebra include(s).
#include "detray/plugins/algebra/array_definitions.hpp"

// Project include(s).
#include "detray/field/constant_magnetic_field.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/simulator.hpp"
#include "tests/common/tools/track_generators.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

using namespace detray;
using transform3 = __plugin::transform3<detray::scalar>;

int main() {

    // Create geometry
    vecmem::host_memory_resource host_mr;
    const auto detector = create_toy_geometry(host_mr);

    // Create B field
    const vector3 b{0, 0, 2 * unit_constants::T};
    constant_magnetic_field<> b_field(b);

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
    auto sim = simulator(n_events, detector, b_field, generator, smearer);

    sim.run();

    return 0;
}