/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/unbounded.hpp"
#include "detray/propagator/detail/jacobian_engine.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/axis_rotation.hpp"
#include "tests/common/tools/intersectors/helix_cylinder_intersector.hpp"
#include "tests/common/tools/intersectors/helix_line_intersector.hpp"
#include "tests/common/tools/intersectors/helix_plane_intersector.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// google-test include(s).
#include <gtest/gtest.h>

using namespace detray;

// Algebra types
using matrix_operator = standard_matrix_operator<scalar>;
using algebra_t = ALGEBRA_PLUGIN<scalar>;
using transform3 = dtransform3D<algebra_t>;
using vector3 = dvector3D<algebra_t>;
using intersection_t =
    intersection2D<surface_descriptor<>, scalar, ALGEBRA_PLUGIN>;

// Mask types to be tested
using rectangle_type = detray::mask<detray::rectangle2D<>>;
using trapezoid_type = detray::mask<detray::trapezoid2D<>>;
using ring_type = detray::mask<detray::ring2D<>>;
using cylinder_type = detray::mask<detray::cylinder2D<>>;
using straw_wire_type = detray::mask<detray::line<false>>;
using cell_wire_type = detray::mask<detray::line<true>>;

// Test class for covariance transport
template <typename T>
class HelixCovarianceTransportValidation : public ::testing::Test {
    public:
    // Environment Setup
    const vector3 B{0.f, 0.f, 1.f * unit<scalar>::T};
    const vector3 y_axis{0.f, 1.f, 0.f};
    const vector3 z_axis{0.f, 0.f, 1.f};
    static constexpr const scalar mask_tolerance{1e-3f};
    static constexpr const scalar tolerance{3e-2f};

    // Test types
    using mask_type = T;
    using local_frame_type = typename mask_type::local_frame_type;
    using helix_intersector_type = helix_intersector<intersection_t, mask_type>;

    // First mask at the origin is always rectangle
    using first_mask_type = rectangle_type;
    using first_jacobian_engine = detail::jacobian_engine<
        first_mask_type::shape::template local_frame_type, scalar,
        ALGEBRA_PLUGIN>;

    // Transform3 type
    using transform3_type = typename local_frame_type::transform3D;
    // Scalar type
    using scalar_type = typename transform3_type::scalar_type;

    // Vector and matrix types
    using jacobian_engine =
        detail::jacobian_engine<mask_type::shape::template local_frame_type,
                                scalar, ALGEBRA_PLUGIN>;
    using free_vector = typename jacobian_engine::free_vector;
    using bound_vector = typename jacobian_engine::bound_vector;
    using bound_matrix = typename jacobian_engine::bound_matrix;
    using free_to_bound_matrix = typename jacobian_engine::free_to_bound_matrix;
    using bound_to_free_matrix = typename jacobian_engine::bound_to_free_matrix;
    using free_matrix = typename jacobian_engine::free_matrix;

    std::tuple<rectangle_type, trapezoid_type, ring_type, cylinder_type,
               straw_wire_type, cell_wire_type>
        masks = std::make_tuple(
            rectangle_type{0u, 50 * unit<scalar>::mm, 50 * unit<scalar>::mm},
            trapezoid_type{0u, 50 * unit<scalar>::mm, 100 * unit<scalar>::mm,
                           50 * unit<scalar>::mm,
                           1.f / (2.f * 50 * unit<scalar>::mm)},
            ring_type{0u, 0.f, 50.f * unit<scalar>::mm},
            cylinder_type{0u, 0.3f * unit<scalar>::mm,
                          -100.f * unit<scalar>::mm, 100.f * unit<scalar>::mm},
            straw_wire_type{0u, 50.f * unit<scalar>::mm,
                            100.f * unit<scalar>::mm},
            cell_wire_type{0u, 50.f * unit<scalar>::mm,
                           100.f * unit<scalar>::mm});

    // Create transform matrices
    std::vector<transform3_type> create_transforms(
        detail::helix<transform3>& reference_helix,
        const std::size_t n_planes) {

        std::vector<transform3> trfs;

        // Step size between two neighbor planes
        const scalar_type S{2.f * constant<scalar>::pi /
                            std::abs(reference_helix._K)};
        const scalar_type step_size = S / static_cast<scalar_type>(n_planes);

        for (std::size_t i = 0u; i < n_planes; i++) {
            const scalar_type s = step_size * scalar_type(i);

            // Translation of the new surface
            vector3 trl = reference_helix(s);

            // Normal vector of the new surface
            vector3 w = reference_helix.dir(s);
            vector3 v = vector::cross(z_axis, w);

            if (i > 0u &&
                (std::is_same_v<local_frame_type, cylindrical2D<algebra_t>> ||
                 std::is_same_v<local_frame_type, line2D<algebra_t>>)) {

                const vector3 r_axis = vector::cross(w, z_axis);

                // Rotate for cylinder and wire
                axis_rotation<transform3> axis_rot(r_axis,
                                                   constant<scalar>::pi / 2.f);

                // Test masks are rotated
                w = axis_rot(w);
                EXPECT_NEAR(getter::norm(r_axis), 1.f, tolerance);

                v = vector::cross(r_axis, w);

                /*
                const vector3 offset{0.f, 0.f, 10.f * unit<scalar>::mm};
                trl = trl + offset;
                */
            } else {

                // @note why does this offset (in y-direction) fail the test???
                // const vector3 offset{0.f, 10.f * unit<scalar>::mm, 10.f *
                // unit<scalar>::mm};
                const vector3 offset{0.f, 0.f, 10.f * unit<scalar>::mm};
                trl = trl + offset;
            }

            // Add transform matrix
            trfs.emplace_back(trl, w, v);
        }

        return trfs;
    }

    // Error propagation
    template <typename departure_mask_type, typename destination_mask_type>
    bound_track_parameters<transform3_type> propagate(
        const bound_track_parameters<transform3>& bound_params,
        const transform3_type& trf_0, const transform3_type& trf_1,
        const departure_mask_type& mask_0, const destination_mask_type& mask_1,
        scalar_type& total_path_length, std::vector<intersection_t>& sfis) {

        using departure_jacobian_engine = detail::jacobian_engine<
            departure_mask_type::shape::template local_frame_type, scalar,
            ALGEBRA_PLUGIN>;

        using destination_jacobian_engine = detail::jacobian_engine<
            destination_mask_type::shape::template local_frame_type, scalar,
            ALGEBRA_PLUGIN>;

        const bound_vector& bound_vec_0 = bound_params.vector();
        const bound_matrix& bound_cov_0 = bound_params.covariance();

        // Free vector at the departure surface
        const auto free_vec_0 = departure_jacobian_engine::bound_to_free_vector(
            trf_0, mask_0, bound_vec_0);

        // Free track at the departure surface
        free_track_parameters<transform3> free_trk_0;
        free_trk_0.set_vector(free_vec_0);

        // Helix from the departure surface
        detail::helix<transform3> hlx(free_trk_0, &B);

        // Bound-to-free jacobian at the departure surface
        const bound_to_free_matrix bound_to_free_jacobi =
            departure_jacobian_engine::bound_to_free_jacobian(trf_0, mask_0,
                                                              bound_vec_0);

        // Get the intersection on the next surface
        const intersection_t is = get_intersection(
            helix_intersector<intersection_t, destination_mask_type>{}(
                hlx, surface_descriptor<>{}, mask_1, trf_1,
                this->mask_tolerance));

        sfis.push_back(is);

        // Helical path length between two surfaces
        const auto path_length = is.path;

        // Add the path length to the total path length
        total_path_length += path_length;

        // Free transport jacobian between two surfaces
        const free_matrix transport_jacobi = hlx.jacobian(path_length);

        // r at the destination surface
        const vector3 r = hlx.pos(path_length);

        // dr/ds, or t at the destination surface
        const vector3 t = hlx.dir(path_length);

        // d^2r/ds^2, or dt/ds at the destination surface
        const vector3 dtds = hlx.qop() * vector::cross(t, B);

        // Free track at the destination surface
        free_track_parameters<transform3> free_trk_1;
        free_trk_1.set_pos(r);
        free_trk_1.set_dir(t);
        free_trk_1.set_qop(free_trk_0.qop());

        // Path correction
        const free_matrix path_correction =
            destination_jacobian_engine::path_correction(r, t, dtds, trf_1);

        // Correction term for the path variation
        const free_matrix correction_term =
            matrix_operator().template identity<e_free_size, e_free_size>() +
            path_correction;

        // Free-to-bound jacobian at the destination surface
        const free_to_bound_matrix free_to_bound_jacobi =
            destination_jacobian_engine::free_to_bound_jacobian(
                trf_1, free_trk_1.vector());

        // Bound vector at the destination surface
        const bound_vector bound_vec_1 =
            destination_jacobian_engine::free_to_bound_vector(
                trf_1, free_trk_1.vector());

        // Full jacobian
        const bound_matrix full_jacobi = free_to_bound_jacobi *
                                         correction_term * transport_jacobi *
                                         bound_to_free_jacobi;

        // Update the covariance at the destination surface
        const bound_matrix bound_cov_1 =
            full_jacobi * bound_cov_0 *
            matrix_operator().transpose(full_jacobi);

        bound_track_parameters<transform3> ret;
        ret.set_vector(bound_vec_1);
        ret.set_covariance(bound_cov_1);

        return ret;
    }

    const intersection_t get_intersection(
        const std::array<intersection_t, 2>& inters) const {
        return inters[0];
    }

    const intersection_t get_intersection(const intersection_t& inters) const {
        return inters;
    }
};

using TestTypes =
    ::testing::Types<rectangle_type, trapezoid_type, ring_type, cylinder_type,
                     straw_wire_type, cell_wire_type>;
TYPED_TEST_SUITE(HelixCovarianceTransportValidation, TestTypes, );

TYPED_TEST(HelixCovarianceTransportValidation, one_loop_test) {

    // @NOTE: The test with high energy (>1 GeV) might fail due
    // to the numerical instability
    free_track_parameters<transform3> free_trk(
        {0.f, 0.f, 0.f}, 0.f, {0.1f * unit<scalar>::GeV, 0.f, 0.f}, -1.f);

    detail::helix<transform3> reference_helix(free_trk, &this->B);

    const std::size_t n_planes = 10u;
    std::vector<transform3> trfs =
        this->create_transforms(reference_helix, n_planes);
    ASSERT_EQ(trfs.size(), 10u);

    // Set the initial bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vec_0 =
        TestFixture::first_jacobian_engine::free_to_bound_vector(
            trfs[0], free_trk.vector());

    // Set the initial bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov_0 =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov_0, e_bound_loc0, e_bound_loc0) = 1.f;
    getter::element(bound_cov_0, e_bound_loc1, e_bound_loc1) = 1.f;
    getter::element(bound_cov_0, e_bound_phi, e_bound_phi) = 1.f;
    // Set theta error to zero, to suppress the loc1 (z) divergence
    getter::element(bound_cov_0, e_bound_theta, e_bound_theta) = 0.f;
    getter::element(bound_cov_0, e_bound_qoverp, e_bound_qoverp) = 1.f;
    getter::element(bound_cov_0, e_bound_time, e_bound_time) = 0.f;

    // Set bound track parameters
    bound_track_parameters<transform3> bound_params;
    bound_params.set_vector(bound_vec_0);
    bound_params.set_covariance(bound_cov_0);

    // Create masks
    const auto first_mask = std::get<rectangle_type>(this->masks);
    const auto test_mask =
        std::get<typename TestFixture::mask_type>(this->masks);

    // Total path length, just for testing purpose
    scalar total_path_length = 0.f;

    // Intersections for testing purporse
    std::vector<intersection_t> sfis;

    // Iterate over the planes until we reach the first plane (one loop)
    for (std::size_t i_p = 0u; i_p < n_planes; i_p++) {

        if (i_p == 0) {
            bound_params =
                this->propagate(bound_params, trfs[i_p], trfs[i_p + 1],
                                first_mask, test_mask, total_path_length, sfis);
        } else if (i_p > 0 && i_p < n_planes - 1) {

            bound_params =
                this->propagate(bound_params, trfs[i_p], trfs[i_p + 1],
                                test_mask, test_mask, total_path_length, sfis);

        } else if (i_p == n_planes - 1) {

            bound_params =
                this->propagate(bound_params, trfs[i_p], trfs[0], test_mask,
                                first_mask, total_path_length, sfis);
        }
    }

    // Check if the total path length is the expected value
    ASSERT_TRUE(total_path_length > 1e-3);
    ASSERT_NEAR(total_path_length,
                2.f * constant<scalar>::pi / std::abs(reference_helix._K),
                this->tolerance);

    ASSERT_EQ(sfis.size(), n_planes);
    for (std::size_t i = 0u; i < n_planes; i++) {
        EXPECT_TRUE(sfis[i].status);
        EXPECT_TRUE(sfis[i].direction);
    }

    // Check if the same vector is obtained after one loop
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec_0, i, 0u),
                    matrix_operator().element(bound_params.vector(), i, 0u),
                    this->tolerance);
    }
    // Check if the same covariance is obtained after one loop
    for (unsigned int i = 0u; i < e_bound_size; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            EXPECT_NEAR(
                matrix_operator().element(bound_cov_0, i, j),
                matrix_operator().element(bound_params.covariance(), i, j),
                this->tolerance);
        }
    }
}
