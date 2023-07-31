/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/utils/consistency_checker.hpp"
#include "tests/common/test_base/fixture_base.hpp"

// System include(s)
#include <iostream>

namespace detray {

/// @brief Test class that runs the consistency check on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class consistency_check : public detray::test::fixture_base<> {
    public:
    using fixture_type = detray::test::fixture_base<>;

    explicit consistency_check(const detector_t &det,
                               const typename detector_t::name_map &names)
        : m_det{det}, m_names{names} {}

    /// Run the consistency check
    void TestBody() override {
        std::cout << "INFO: Running consistency check on: " << m_names.at(0)
                  << std::endl;
        ASSERT_TRUE(detail::check_consistency(m_det));
    }

    private:
    /// The detector to be checked
    const detector_t &m_det;
    /// Volume names
    const typename detector_t::name_map &m_names;
};

}  // namespace detray
