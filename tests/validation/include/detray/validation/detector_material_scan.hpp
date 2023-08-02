/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/io/common/detail/file_handle.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/types.hpp"
#include "tests/common/test_base/fixture_base.hpp"
#include "tests/common/tools/particle_gun.hpp"

// System include(s)
#include <ios>
#include <iostream>
#include <string>

namespace detray {

/// @brief Test class that runs the material ray scan on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class material_scan : public test::fixture_base<> {

    using transform3_t = typename detector_t::transform3;
    using scalar_t = typename detector_t::scalar_type;
    using ray_t = detail::ray<transform3_t>;

    public:
    using fixture_type = test::fixture_base<>;

    struct config : public fixture_type::configuration {
        using trk_gen_config_t =
            typename uniform_track_generator<ray_t>::configuration;

        std::string m_name{"material_scan"};
        trk_gen_config_t m_trk_gen_cfg{};

        /// Getters
        /// @{
        const std::string &name() const { return m_name; }
        trk_gen_config_t &track_generator() { return m_trk_gen_cfg; }
        const trk_gen_config_t &track_generator() const {
            return m_trk_gen_cfg;
        }
        /// @}

        /// Setters
        /// @{
        config &name(const std::string n) {
            m_name = n;
            return *this;
        }
        /// @}
    };

    template <typename config_t>
    explicit material_scan(const detector_t &det,
                           const typename detector_t::name_map &names,
                           const config_t &cfg = {})
        : m_det{det}, m_names{names} {
        m_cfg.name(cfg.name());
        m_cfg.track_generator() = cfg.track_generator();
    }

    /// Run the ray scan
    void TestBody() override {

        std::size_t n_tracks{0u};
        auto ray_generator =
            uniform_track_generator<ray_t>(m_cfg.track_generator());

        // Csv output file
        std::string file_name{m_cfg.name() + "_" + m_names.at(0)};
        detray::io::detail::file_handle outfile{
            file_name, ".csv",
            std::ios::out | std::ios::binary | std::ios::trunc};
        *outfile << "eta,phi,mat_sX0,mat_sL0,mat_tX0,mat_tL0" << std::endl;

        std::cout << "INFO: Running material scan on: " << m_names.at(0)
                  << "\n(" << ray_generator.size() << " rays) ...\n"
                  << std::endl;

        scalar_t eta{}, phi{}, mat_sX0{}, mat_sL0{}, mat_tX0{}, mat_tL0{};
        for (const auto ray : ray_generator) {

            // Record all intersections and surfaces along the ray
            const auto intersection_record =
                particle_gun::shoot_particle(m_det, ray);

            if (intersection_record.empty()) {
                std::cout << "ERROR: Intersection trace empty for ray "
                          << n_tracks << "/" << ray_generator.size() << ": "
                          << ray << std::endl;
                break;
            }

            eta = getter::eta(ray.dir());
            phi = getter::phi(ray.dir());
            // Total accumulated path length in X0 and L0, respectively
            mat_sX0 = 0.f;
            mat_sL0 = 0.f;
            // Total accumulated material thickness in X0 and L0, respectively
            mat_tX0 = 0.f;
            mat_tL0 = 0.f;

            // Record material for this ray
            for (const auto &record : intersection_record) {

                const auto sf = surface{m_det, record.second.sf_desc};

                const auto [seg, t, mx0, ml0] =
                    sf.template visit_material<get_material_params>(
                        record.second.local, record.second.cos_incidence_angle);

                if (mx0 > 0.f) {
                    mat_sX0 += seg / mx0;
                    mat_tX0 += t / mx0;
                } else {
                    std::cout << "WARNING: Encountered invalid X_0: " << mx0
                              << "\nOn surface: " << sf << std::endl;
                }
                if (ml0 > 0.f) {
                    mat_sL0 += seg / ml0;
                    mat_tL0 += t / ml0;
                } else {
                    std::cout << "WARNING: Encountered invalid L_0: " << ml0
                              << "\nOn surface: " << sf << std::endl;
                }
            }

            if (mat_sX0 == 0.f or mat_sL0 == 0.f or mat_tX0 == 0.f or
                mat_tL0 == 0.f) {
                std::cout << "WARNING: No material recorded for ray "
                          << n_tracks << "/" << ray_generator.size() << ": "
                          << ray << std::endl;
            }

            *outfile << eta << "," << phi << "," << mat_sX0 << "," << mat_sL0
                     << "," << mat_tX0 << "," << mat_tL0 << std::endl;

            ++n_tracks;
        }
    }

    private:
    /// @brief Functor to retrieve the material parameters of a given
    /// intersection
    struct get_material_params {

        /// Access to material slabs or rods in a homogeneous material
        /// description
        template <class material_coll_t, class point_t,
                  std::enable_if_t<not detail::is_grid_v<
                                       typename material_coll_t::value_type>,
                                   bool> = true>
        inline constexpr decltype(auto) get_material(
            const material_coll_t &material_coll, const dindex idx,
            const point_t &) const noexcept {
            return material_coll[idx];
        }

        /// Access to material slabs in a material map
        template <class material_coll_t, class point_t,
                  std::enable_if_t<
                      detail::is_grid_v<typename material_coll_t::value_type>,
                      bool> = true>
        inline constexpr decltype(auto) get_material(
            const material_coll_t &material_coll, const dindex idx,
            const point_t &loc_point) const noexcept {

            // Find the material slab (only one entry per bin)
            return *(material_coll[idx].search(loc_point));
        }

        template <typename mat_group_t, typename index_t, typename point_t>
        inline auto operator()(const mat_group_t &mat_group,
                               const index_t &index, const point_t &loc,
                               const scalar_t cos_inc_angle) const {

            const auto slab = get_material(mat_group, index, loc);

            const scalar_t seg{slab.path_segment(cos_inc_angle, loc[0])};
            const scalar_t t{slab.thickness()};
            const scalar_t mat_X0{slab.get_material().X0()};
            const scalar_t mat_L0{slab.get_material().L0()};

            return std::tuple(seg, t, mat_X0, mat_L0);
        }
    };

    /// The configuration of this test
    config m_cfg;
    /// The detector to be checked
    const detector_t &m_det;
    /// Volume names
    const typename detector_t::name_map &m_names;
};

}  // namespace detray
