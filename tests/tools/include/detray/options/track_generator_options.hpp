/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"

// Detray test include(s)
#include "detray/options/options_handling.hpp"
#include "detray/test/utils/simulation/event_generator/random_track_generator_config.hpp"
#include "detray/test/utils/simulation/event_generator/uniform_track_generator_config.hpp"

// Boost
#include "detray/options/boost_program_options.hpp"

// System include(s)
#include <stdexcept>
#include <vector>

namespace detray::options {

namespace detail {

/// Add options for detray event generation
template <concepts::scalar scalar_t>
void add_uniform_track_gen_options(
    boost::program_options::options_description &desc,
    const uniform_track_generator_config<scalar_t> &cfg) {

    desc.add_options()(
        "phi_steps",
        boost::program_options::value<std::size_t>()->default_value(
            cfg.phi_steps()),
        "No. phi steps for particle gun")(
        "eta_steps",
        boost::program_options::value<std::size_t>()->default_value(
            cfg.eta_steps()),
        "No. eta steps for particle gun")(
        "eta_range",
        boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
        "Min, Max range of eta values for particle gun")(
        "randomize_charge", "Randomly flip charge sign per track")(
        "origin",
        boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
        "Coordintates for particle gun origin position [mm]")(
        "p_tot",
        boost::program_options::value<scalar_t>()->default_value(
            static_cast<scalar_t>(cfg.m_p_mag) / unit<scalar_t>::GeV),
        "Total momentum of the test particle [GeV]")(
        "p_T",
        boost::program_options::value<scalar_t>()->default_value(
            static_cast<scalar_t>(cfg.m_p_mag) / unit<scalar_t>::GeV),
        "Transverse momentum of the test particle [GeV]");
}

/// Add options for detray event generation
template <concepts::scalar scalar_t>
void configure_uniform_track_gen_options(
    const boost::program_options::variables_map &vm,
    uniform_track_generator_config<scalar_t> &cfg) {

    cfg.phi_steps(vm["phi_steps"].as<std::size_t>());
    cfg.eta_steps(vm["eta_steps"].as<std::size_t>());
    cfg.randomize_charge(vm.count("randomize_charge"));

    if (vm.count("eta_range")) {
        const auto eta_range = vm["eta_range"].as<std::vector<scalar_t>>();
        if (eta_range.size() == 2u) {
            cfg.eta_range(eta_range[0], eta_range[1]);
        } else {
            throw std::invalid_argument("Eta range needs two arguments");
        }
    }
    if (vm.count("origin")) {
        const auto origin = vm["origin"].as<std::vector<scalar_t>>();
        if (origin.size() == 3u) {
            cfg.origin(origin[0] * unit<scalar_t>::mm,
                       origin[1] * unit<scalar_t>::mm,
                       origin[2] * unit<scalar_t>::mm);
        } else {
            throw std::invalid_argument(
                "Particle gun origin needs three arguments");
        }
    }
    if (!vm["p_T"].defaulted() && !vm["p_tot"].defaulted()) {
        throw std::invalid_argument(
            "Transverse and total momentum cannot be specified at the same "
            "time");
    }
    if (!vm["p_T"].defaulted()) {
        cfg.p_T(vm["p_T"].as<scalar_t>() * unit<scalar_t>::GeV);
    } else {
        cfg.p_tot(vm["p_tot"].as<scalar_t>() * unit<scalar_t>::GeV);
    }
}

/// Add options for detray event generation
template <concepts::scalar scalar_t>
void add_rnd_track_gen_options(
    boost::program_options::options_description &desc,
    const random_track_generator_config<scalar_t> &cfg) {

    desc.add_options()(
        "n_tracks",
        boost::program_options::value<std::size_t>()->default_value(
            cfg.n_tracks()),
        "No. of tracks for particle gun")(
        "theta_range",
        boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
        "Min, Max range of theta values for particle gun")(
        "eta_range",
        boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
        "Min, Max range of eta values for particle gun")(
        "randomize_charge", "Randomly flip charge sign per track")(
        "origin",
        boost::program_options::value<std::vector<scalar_t>>()->multitoken(),
        "Coordintates for particle gun origin position")(
        "p_tot",
        boost::program_options::value<scalar_t>()->default_value(
            static_cast<scalar_t>(cfg.mom_range()[0]) / unit<scalar_t>::GeV),
        "Total momentum of the test particle [GeV]")(
        "p_T",
        boost::program_options::value<scalar_t>()->default_value(
            static_cast<scalar_t>(cfg.mom_range()[0]) / unit<scalar_t>::GeV),
        "Transverse momentum of the test particle [GeV]");
}

/// Add options for detray event generation
template <concepts::scalar scalar_t>
void configure_rnd_track_gen_options(
    const boost::program_options::variables_map &vm,
    random_track_generator_config<scalar_t> &cfg) {

    cfg.n_tracks(vm["n_tracks"].as<std::size_t>());
    cfg.randomize_charge(vm.count("randomize_charge"));

    if (vm.count("eta_range") && vm.count("theta_range")) {
        throw std::invalid_argument(
            "Eta range and theta range cannot be specified at the same time");
    } else if (vm.count("eta_range")) {
        const auto eta_range = vm["eta_range"].as<std::vector<scalar_t>>();
        if (eta_range.size() == 2u) {
            scalar_t min_theta{2.f * std::atan(std::exp(-eta_range[0]))};
            scalar_t max_theta{2.f * std::atan(std::exp(-eta_range[1]))};

            // Wrap around
            if (min_theta > max_theta) {
                scalar_t tmp{min_theta};
                min_theta = max_theta;
                max_theta = tmp;
            }

            cfg.theta_range(min_theta, max_theta);
        } else {
            throw std::invalid_argument("Eta range needs two arguments");
        }
    } else if (vm.count("theta_range")) {
        const auto theta_range = vm["theta_range"].as<std::vector<scalar_t>>();
        if (theta_range.size() == 2u) {
            cfg.theta_range(theta_range[0], theta_range[1]);
        } else {
            throw std::invalid_argument("Theta range needs two arguments");
        }
    }
    if (vm.count("origin")) {
        const auto origin = vm["origin"].as<std::vector<scalar_t>>();
        if (origin.size() == 3u) {
            cfg.origin(origin[0] * unit<scalar_t>::mm,
                       origin[1] * unit<scalar_t>::mm,
                       origin[2] * unit<scalar_t>::mm);
        } else {
            throw std::invalid_argument(
                "Particle gun origin needs three coordinates");
        }
    }
    if (!vm["p_T"].defaulted() && !vm["p_tot"].defaulted()) {
        throw std::invalid_argument(
            "Transverse and total momentum cannot be specified at the same "
            "time");
    }
    if (!vm["p_T"].defaulted()) {
        cfg.p_T(vm["p_T"].as<scalar_t>() * unit<scalar_t>::GeV);
    } else {
        cfg.p_tot(vm["p_tot"].as<scalar_t>() * unit<scalar_t>::GeV);
    }
}

}  // namespace detail

/// Add options for the uniform track generator
/// @{
template <>
void add_options<uniform_track_generator_config<float>>(
    boost::program_options::options_description &desc,
    const uniform_track_generator_config<float> &cfg) {
    detail::add_uniform_track_gen_options(desc, cfg);
}

template <>
void add_options<uniform_track_generator_config<double>>(
    boost::program_options::options_description &desc,
    const uniform_track_generator_config<double> &cfg) {
    detail::add_uniform_track_gen_options(desc, cfg);
}
/// @}

/// Configure the detray uniform track generator
/// @{
template <>
void configure_options<uniform_track_generator_config<float>>(
    const boost::program_options::variables_map &vm,
    uniform_track_generator_config<float> &cfg) {
    detail::configure_uniform_track_gen_options(vm, cfg);
}

template <>
void configure_options<uniform_track_generator_config<double>>(
    const boost::program_options::variables_map &vm,
    uniform_track_generator_config<double> &cfg) {
    detail::configure_uniform_track_gen_options(vm, cfg);
}
/// @}

/// Add options for the random track generator
/// @{
template <>
void add_options<random_track_generator_config<float>>(
    boost::program_options::options_description &desc,
    const random_track_generator_config<float> &cfg) {
    detail::add_rnd_track_gen_options(desc, cfg);
}

template <>
void add_options<random_track_generator_config<double>>(
    boost::program_options::options_description &desc,
    const random_track_generator_config<double> &cfg) {
    detail::add_rnd_track_gen_options(desc, cfg);
}
/// @}

/// Configure the detray random track generator
/// @{
template <>
void configure_options<random_track_generator_config<float>>(
    const boost::program_options::variables_map &vm,
    random_track_generator_config<float> &cfg) {
    detail::configure_rnd_track_gen_options(vm, cfg);
}

template <>
void configure_options<random_track_generator_config<double>>(
    const boost::program_options::variables_map &vm,
    random_track_generator_config<double> &cfg) {
    detail::configure_rnd_track_gen_options(vm, cfg);
}
/// @}

}  // namespace detray::options
