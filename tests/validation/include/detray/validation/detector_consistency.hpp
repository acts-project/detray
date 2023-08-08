/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/volume_graph.hpp"
#include "detray/utils/consistency_checker.hpp"
#include "tests/common/test_base/fixture_base.hpp"
#include "tests/common/tools/hash_tree.hpp"

// System include(s)
#include <iostream>
#include <string>

namespace detray {

/// @brief Test class that runs the consistency check on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class consistency_check : public detray::test::fixture_base<> {
    public:
    using fixture_type = detray::test::fixture_base<>;

    struct config : public fixture_type::configuration {
        std::string m_name{"detector_consistency"};
        bool m_write_graph{false};

        /// Getters
        /// @{
        const std::string &name() const { return m_name; }
        bool write_graph() const { return m_write_graph; }
        /// @}

        /// Setters
        /// @{
        config &name(const std::string n) {
            m_name = n;
            return *this;
        }
        config &write_graph(const bool do_write) {
            m_write_graph = do_write;
            return *this;
        }
        /// @}
    };

    template <typename config_t>
    explicit consistency_check(const detector_t &det,
                               const typename detector_t::name_map &names,
                               const config_t &cfg = {})
        : m_cfg{cfg}, m_det{det}, m_names{names} {}

    /// Run the consistency check
    void TestBody() override {
        std::cout << "INFO: Running consistency check on: " << m_names.at(0)
                  << std::endl
                  << std::endl;

        // Build the graph
        volume_graph graph(m_det);

        ASSERT_TRUE(detail::check_consistency(m_det)) << graph.to_string();

        if (m_cfg.write_graph()) {
            std::cout << graph.to_string() << std::endl;
        }

        // Not currently supported
        if (false) {
            constexpr std::size_t root_hash{3244u};  // < toy detector only

            const auto &adj_mat = graph.adjacency_matrix();

            auto geo_checker = hash_tree(adj_mat);

            EXPECT_EQ(geo_checker.root(), root_hash) << graph.to_string();
        }
    }

    private:
    /// The configuration of this test
    const config m_cfg;
    /// The detector to be checked
    const detector_t &m_det;
    /// Volume names
    const typename detector_t::name_map &m_names;
};

}  // namespace detray
