// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <ostream>
#include <string>

namespace detray::benchmarks {

/// Base type for detray benchmarks with google benchmark
struct benchmark_base {
    /// Local configuration type
    struct configuration {
        /// Size of data sample to be used in benchmark
        int m_samples{100};
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
        constexpr int n_samples() const { return m_samples; }
        constexpr bool do_warmup() const { return m_warmup; }
        constexpr int n_warmup() const { return m_n_warmup; }
        /// @}

        private:
        /// Print the benchmark setup
        friend std::ostream& operator<<(std::ostream& os,
                                        const configuration& cfg) {
            os << " -> running:\t " << cfg.n_samples() << " samples"
               << std::endl;
            if (cfg.do_warmup()) {
                os << " -> warmup: \t " << cfg.n_warmup() << " samples"
                   << std::endl;
            }
            os << std::endl;
            return os;
        }
    };

    /// Default construction
    benchmark_base() = default;

    /// Default destructor
    virtual ~benchmark_base() = default;
};

}  // namespace detray::benchmarks
