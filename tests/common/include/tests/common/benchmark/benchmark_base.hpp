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
        std::size_t m_samples{100u};
        /// Run a number of operations before the benchmark
        bool m_warmup = true;
        // Sleep after building data sample
        bool m_sleep = false;
        // Size of data in warm-up round
        std::size_t m_n_warmup{
            static_cast<std::size_t>(0.1 * static_cast<double>(m_samples))};
        // Size of data in warm-up round
        std::size_t m_n_sleep{1u};

        /// Setters
        /// @{
        configuration& n_samples(std::size_t n) {
            m_samples = n;
            return *this;
        }
        configuration& do_warmup(bool b) {
            m_warmup = b;
            return *this;
        }
        configuration& n_warmup(std::size_t n) {
            m_n_warmup = n;
            m_warmup = true;
            return *this;
        }
        configuration& do_sleep(bool b) {
            m_sleep = b;
            return *this;
        }
        configuration& n_sleep(std::size_t n) {
            m_n_sleep = n;
            m_sleep = true;
            return *this;
        }
        /// @}

        /// Getters
        /// @{
        std::size_t n_samples() const { return m_samples; }
        constexpr bool do_warmup() const { return m_warmup; }
        constexpr std::size_t n_warmup() const { return m_n_warmup; }
        constexpr bool do_sleep() const { return m_sleep; }
        constexpr std::size_t n_sleep() const { return m_n_sleep; }
        /// @}

        /// Print configuration
        friend std::ostream& operator<<(std::ostream& os,
                                        const configuration& c);
    };

    /// Default construction
    benchmark_base() = default;

    /// Default destructor
    virtual ~benchmark_base() = default;

    /// @returns the benchmark name
    virtual std::string name() const = 0;
};

std::ostream& operator<<(std::ostream& os,
                         const benchmark_base::configuration& cfg) {
    os << " -> running:\t " << cfg.n_samples() << " samples" << std::endl;
    if (cfg.do_warmup()) {
        os << " -> warmup: \t " << cfg.n_warmup() << " samples" << std::endl;
    }
    if (cfg.do_sleep()) {
        os << " -> cool down:\t " << cfg.n_sleep() << "s" << std::endl;
    }
    os << std::endl;
    return os;
}

}  // namespace detray