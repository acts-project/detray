/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)#include "detray/definitions/units.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "tests/common/test_base/fixture_base.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <cassert>

namespace detray {

// @todo: make bfield configurable when its available
/*template <typename bfield_bknd_t>*/
struct toy_detector_fixture : public test::fixture_base {
    /// Useful typedefs
    /// @{
    using scalar_type = detray::scalar;
    using point2_type = __plugin::point2<scalar>;
    using point3_type = __plugin::point3<scalar>;
    using vector3_type = __plugin::vector3<scalar>;
    using transform3_type = __plugin::transform3<detray::scalar>;
    using matrix_operator = typename transform3_type::matrix_actor;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;
    /// @}

    /// Detector dependent types
    /// @{
    using detector_type = detector<detector_registry::toy_detector>;
    using volume_type = typename detector_type::volume_type;
    using surface_type = typename detector_type::surface_type;
    using bfield_type = typename detector_type::bfield_type;
    /// @}

    /// Prefix for the benchmark name
    inline static const std::string s_name{"toy_detector"};

    /// Local configuration type
    struct configuration : test::fixture_base::configuration {
        /// Number of pixel barrel layers
        unsigned int m_n_brl_layers{4u};
        /// Number of pixel endcap layers on either side
        unsigned int m_n_edc_layers{3u};

        /// Default construciton
        configuration() = default;

        /// Construct from a base configuration
        configuration(const detray::test::fixture_base::configuration &cfg)
            : detray::test::fixture_base::configuration(cfg) {}

        configuration &n_barrel_layers(unsigned int n) {
            m_n_brl_layers = n;
            return *this;
        }

        configuration &n_endcap_layers(unsigned int n) {
            m_n_edc_layers = n;
            return *this;
        }

        unsigned int n_barrel_layers() const { return m_n_brl_layers; }
        unsigned int n_endcap_layers() const { return m_n_edc_layers; }
    };

    /// The detector configuration
    configuration m_cfg{};

    /// Default construction
    toy_detector_fixture() = default;

    /// Construct from an externally provided configuration @param cfg
    toy_detector_fixture(configuration cfg) : m_cfg{cfg} {}

    /// @return the benchmark configuration
    configuration &config() { return m_cfg; }

    /// Give better names to tests and benchmarks
    std::string name() const { return toy_detector_fixture::s_name; }

    /// Prepare data and run benchmark loop
    detector_type build_detector(vecmem::memory_resource *mr) const {

        // Create detector
        return create_toy_geometry(*mr, m_cfg.n_barrel_layers(),
                                   m_cfg.n_endcap_layers());
    }
};

}  // namespace detray