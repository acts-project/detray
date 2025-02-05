// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Boost
#include "detray/options/boost_program_options.hpp"

// System include(s)
#include <iostream>

namespace detray::options {

/// Add options to the boost options description according to the type T
template <typename T>
void add_options(boost::program_options::options_description &,
                 const T &) { /* Do nothing */
}

/// Fill the configuration type T from the boost variable map
template <typename T>
void configure_options(const boost::program_options::variables_map &,
                       T &) { /* Do nothing */
}

/// Print the configuration
template <typename T>
void print_options(T &cfg) {
    std::cout << cfg << std::endl;
}

}  // namespace detray::options
