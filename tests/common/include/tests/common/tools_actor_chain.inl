/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <sstream>
#include <string>

#include "detray/definitions/units.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/base_actor.hpp"

namespace {

using namespace detray;
using namespace __plugin;

/// Actor that prints its ID and a number
struct print_actor : detray::actor {

    struct print_actor_state {
        std::stringstream stream{};

        std::string to_string() const { return stream.str(); }
    };

    // Broadcast state type to actor chain
    using state_type = print_actor_state;

    /// Print the ID of the type and the id in its state
    template <typename propagator_state_t>
    void operator()(state_type &printer_state,
                    const propagator_state_t & /*p_state*/) const {
        printer_state.stream << "[print actor]:";
    }

    template <typename subj_state_t, typename propagator_state_t>
    void operator()(state_type &printer_state,
                    const subj_state_t &subject_state,
                    const propagator_state_t & /*p_state*/) const {
        printer_state.stream << "[print actor obs "
                             << subject_state.some_state.back() << "]:";
    }
};

/// Example actor that
template <template <typename...> class vector_t>
struct example_actor : detray::actor {

    /// actor state
    struct example_actor_state {

        example_actor_state() {}

        // Keep dynamic data per propagation stream
        vector_t<float> some_state;
    };

    // Broadcast state type to actor chain
    using state_type = example_actor_state;

    /// Actor implementation
    template <typename propagator_state_t>
    void operator()(state_type &example_state,
                    const propagator_state_t & /*p_state*/) const {
        example_state.some_state.push_back(example_state.some_state.size());
    }

    /// Actor implementation when observing another actor
    template <typename subj_state_t, typename propagator_state_t>
    void operator()(state_type &example_state,
                    const subj_state_t &subject_state,
                    const propagator_state_t & /*p_state*/) const {
        example_state.some_state.push_back(subject_state.some_state.size() /
                                           10.);
    }
};

// Implements example_actor with two print observers
using example_actor_t = example_actor<std::vector>;
using composite1 =
    composite_actor<dtuple, example_actor_t, print_actor, print_actor>;
// Implements example_actor with one print observer
using composite2 = composite_actor<dtuple, example_actor_t, print_actor>;
// Implements example_actor through composite2 and has composite1 as observer
using composite3 = composite_actor<dtuple, composite2, composite1>;
// Implements example_actor through composite2<-composite3 with composite1 obs.
using composite4 = composite_actor<dtuple, composite3, composite1>;

}  // anonymous namespace

// Test the actor chain on some dummy actor types
TEST(ALGEBRA_PLUGIN, actor_chain) {

    // The actor states (can be reused between actors, e.g. for the printer or
    // empty states)
    example_actor_t::state_type example_state{};
    print_actor::state_type printer_state{};

    // Aggregate actor states to be able to pass them through the chain
    auto actor_states = std::tie(example_state, printer_state);

    // Propagator state
    struct empty_prop_state {};
    empty_prop_state prop_state{};

    // Chain of actors
    using actor_chain_t = actor_chain<dtuple, example_actor_t, composite1,
                                      composite2, composite3, composite4>;
    // Run
    actor_chain_t run_actors{};
    run_actors(actor_states, prop_state);

    ASSERT_TRUE(printer_state.to_string().compare(
                    "[print actor obs 1]:[print actor obs 1]:[print actor obs "
                    "2]:[print actor obs 0.4]:[print actor obs 0.4]:[print "
                    "actor obs 0.6]:[print actor obs 0.6]:") == 0)
        << "Printer call chain: " << printer_state.to_string() << std::endl;
}
