/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/benchmarks/benchmark_base.hpp"
#include "detray/propagator/propagation_config.hpp"

// System include(s)
#include <string>

namespace detray {

/// Configuration for propagation benchmarks
struct propagation_benchmark_config {
    /// Prefix for the benchmark name
    std::string m_name{"BM_PROPAGATION"};
    /// Benchmark configuration
    benchmark_base::configuration m_benchmark{};
    /// Propagation configuration
    propagation::config m_propagation{};

    /// Default construciton
    propagation_benchmark_config() = default;

    /// Construct from a base configuration
    propagation_benchmark_config(
        const detray::benchmark_base::configuration& bench_cfg)
        : m_benchmark(bench_cfg) {}

    /// Getters
    /// @{
    const std::string& name() const { return m_name; }
    const propagation::config& propagation() const { return m_propagation; }
    propagation::config& propagation() { return m_propagation; }
    const benchmark_base::configuration& benchmark() const {
        return m_benchmark;
    }
    benchmark_base::configuration& benchmark() { return m_benchmark; }
    /// @}

    /// Setters
    /// @{
    propagation_benchmark_config& name(std::string& n) {
        m_name = n;
        return *this;
    }
    /// @}
};

}  // namespace detray
