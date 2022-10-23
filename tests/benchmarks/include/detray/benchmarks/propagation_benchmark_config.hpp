/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/benchmarks/benchmark_base.hpp"
#include "detray/propagator/propagation_config.hpp"

// System include(s)
#include <string>
#include <string_view>

namespace detray::benchmarks {

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
    explicit propagation_benchmark_config(
        const benchmark_base::configuration& bench_cfg)
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
    propagation_benchmark_config& name(const std::string_view n) {
        m_name = n;
        return *this;
    }
    /// @}
};

}  // namespace detray::benchmarks
