/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <ostream>
#include <string>

namespace detray {

/// Base type for linear algebra benchmarks with google benchmark
struct benchmark_base {
    /// Local configuration type
    struct configuration {
        /// Size of data sample to be used in benchmark
        int m_samples{100u};
        /// Run a number of operations before the benchmark
        bool m_warmup = true;
        // Size of data in warm-up round
        int m_n_warmup{static_cast<int>(0.1 * static_cast<double>(m_samples))};

        /// Setters
        /// @{
        configuration& n_samples(int n) {
            m_samples = n;
            return *this;
        }
        configuration& do_warmup(bool b) {
            m_warmup = b;
            return *this;
        }
        configuration& n_warmup(int n) {
            m_n_warmup = n;
            m_warmup = true;
            return *this;
        }
        /// @}

        /// Getters
        /// @{
        int n_samples() const { return m_samples; }
        constexpr bool do_warmup() const { return m_warmup; }
        constexpr int n_warmup() const { return m_n_warmup; }
        /// @}

        /// Print configuration
        friend std::ostream& operator<<(std::ostream& os,
                                        const configuration& c);
    };

    /// Default construction
    benchmark_base() = default;

    /// Default destructor
    virtual ~benchmark_base() = default;
};

inline std::ostream& operator<<(std::ostream& os,
                                const benchmark_base::configuration& cfg) {
    os << " -> running:\t " << cfg.n_samples() << " samples" << std::endl;
    if (cfg.do_warmup()) {
        os << " -> warmup: \t " << cfg.n_warmup() << " samples" << std::endl;
    }
    os << std::endl;
    return os;
}

}  // namespace detray
