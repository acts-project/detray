/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s).
#include "sf_finders_grid_cuda_kernel.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <limits>

using namespace detray;

namespace {

/// Test single entry bin content element by element
template <typename content_t>
void test_content(const content_t& bin_content, const content_t& expected) {
    unsigned int i{0u};
    for (const auto& elem : bin_content) {
        ASSERT_NEAR(elem, expected[i++], tol);
    }
}

/// Test collection entries in a bin element by element
template <typename content_t, typename expected_content_t>
void test_entry_collection(const content_t& bin_content,
                           const expected_content_t& expected) {

    // Running indices for the points in the collection and the elements of a
    // single point3
    unsigned int i = 0u, j = 0u;
    for (const auto& entry : bin_content) {
        const auto& expt_entry = expected[i++];
        j = 0u;
        for (const auto& elem : entry) {
            ASSERT_NEAR(elem, expt_entry[j++], tol);
        }
    }
}

// Create some bin data for non-owning grid
template <class populator_t, typename entry_t>
struct bin_content_sequence {

    entry_t entry{0};

    auto operator()() {
        entry += entry_t{1};
        return populator_t::init(entry);
    }
};

}  // anonymous namespace

TEST(grids_cuda, grid3_replace_populator) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // Build multi-axis
    using axes_t = cartesian_3D<host_container_types>;

    typename axes_t::boundary_storage_type axis_data(&mng_mr);
    typename axes_t::edges_storage_type bin_edges(&mng_mr);

    axis_data.reserve(3);
    axis_data.insert(
        axis_data.begin(),
        {dindex_range{0u, 3u}, dindex_range{2u, 6u}, dindex_range{4u, 10u}});
    bin_edges.reserve(6);
    bin_edges.insert(bin_edges.begin(), {-1.f, 4.f, 0.f, 6.f, -5.f, 5.f});

    axes_t axes(std::move(axis_data), std::move(bin_edges));

    // build host grid
    host_grid3_replace::bin_storage_type bin_data(&mng_mr);
    bin_data.resize(
        3u * 6u * 10u,
        populator<host_grid3_replace::populator_impl>::init<point3>());

    host_grid3_replace g3(std::move(bin_data), std::move(axes));

    const auto& axis_x = g3.template get_axis<n_axis::label::e_x>();
    const auto& axis_y = g3.template get_axis<n_axis::label::e_y>();
    const auto& axis_z = g3.template get_axis<n_axis::label::e_z>();

    // pre-check
    for (unsigned int i_x = 0u; i_x < axis_x.nbins(); i_x++) {
        for (unsigned int i_y = 0u; i_y < axis_y.nbins(); i_y++) {
            for (unsigned int i_z = 0u; i_z < axis_z.nbins(); i_z++) {
                const auto& bin = g3.bin({i_x, i_y, i_z});
                auto invalid_bin = populator<
                    host_grid3_replace::populator_impl>::init<point3>();
                test_content(bin[0], invalid_bin.content());
            }
        }
    }
    // run device side test, which populates the grid
    grid_replace_test(get_data(g3), axis_x.nbins(), axis_y.nbins(),
                      axis_z.nbins());
    // post-check
    for (unsigned int i_x = 0u; i_x < axis_x.nbins(); i_x++) {
        for (unsigned int i_y = 0u; i_y < axis_y.nbins(); i_y++) {
            for (unsigned int i_z = 0u; i_z < axis_z.nbins(); i_z++) {
                const dindex gbin_idx =
                    g3.serialize(detray::n_axis::multi_bin<3>{i_x, i_y, i_z});

                const auto& bin = g3.bin(gbin_idx);

                const detray::scalar gbin_idx_f{
                    static_cast<detray::scalar>(gbin_idx)};
                const point3 tp{axis_x.min() + gbin_idx_f * axis_x.bin_width(),
                                axis_y.min() + gbin_idx_f * axis_y.bin_width(),
                                axis_z.min() + gbin_idx_f * axis_z.bin_width()};

                test_content(bin[0], tp);
            }
        }
    }
}

TEST(grids_cuda, grid2_replace_populator_ci) {

    vecmem::cuda::managed_memory_resource mng_mr;

    // Build multi-axis
    using axes_t = polar_ir<host_container_types>;

    typename axes_t::boundary_storage_type axis_data(&mng_mr);
    typename axes_t::edges_storage_type bin_edges(&mng_mr);

    axis_data.reserve(2u);
    axis_data.insert(axis_data.end(),
                     {dindex_range{0u, 4u}, dindex_range{5u, 2u}});
    bin_edges.reserve(7u);
    bin_edges.insert(bin_edges.end(), {1.f, 3.f, 9.f, 27.f, 81.f, -2.f, 4.f});

    axes_t axes(std::move(axis_data), std::move(bin_edges));

    // build host grid
    host_grid2_replace_ci::bin_storage_type bin_data(&mng_mr);
    bin_data.resize(
        2u * 4u,
        populator<host_grid2_replace_ci::populator_impl>::init<point3>());

    host_grid2_replace_ci g2(std::move(bin_data), std::move(axes));

    const auto& axis_r = g2.template get_axis<n_axis::label::e_r>();
    const auto& axis_phi = g2.template get_axis<n_axis::label::e_phi>();

    // pre-check
    for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
        for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
            const auto& bin = g2.bin({i_r, i_phi});
            auto invalid_bin = populator<
                host_grid2_replace_ci::populator_impl>::init<point3>();

            test_content(bin[0], invalid_bin.content());
        }
    }

    // run device side test, which populates the grid
    grid_replace_ci_test(get_data(g2), axis_r.nbins(), axis_phi.nbins());

    // post-check
    for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
        for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
            const dindex gbin_idx =
                g2.serialize(detray::n_axis::multi_bin<2>{i_r, i_phi});
            const auto& bin = g2.bin(gbin_idx);

            const detray::scalar gbin_idx_f{
                static_cast<detray::scalar>(gbin_idx)};
            const point3 tp{axis_r.min() + gbin_idx_f * axis_r.bin_width(i_r),
                            axis_phi.min() + gbin_idx_f * axis_phi.bin_width(),
                            0.5f};

            test_content(bin[0], tp);
        }
    }
}

/*TEST(grids_cuda, grid2_complete_populator) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // Build multi-axis
    using axes_t = polar<host_container_types>;

    typename axes_t::boundary_storage_type axis_data(&mng_mr);
    typename axes_t::edges_storage_type bin_edges(&mng_mr);

    axis_data.reserve(2u);
    axis_data.insert(axis_data.begin(),
                      {dindex_range{0u, 3u}, dindex_range{2u, 7u}});
    bin_edges.reserve(4);
    bin_edges.insert(bin_edges.begin(), {0.f, 3.f, -1.f, 6.f});

    axes_t axes(std::move(axis_data), std::move(bin_edges));

    // build host grid
    const point3 first_tp{3.f, 3.f, 3.f};

    host_grid2_complete::bin_storage_type bin_data(&mng_mr);
    bin_data.resize(
        3u * 7u,
populator<host_grid2_complete::populator_impl>::init<point3>(first_tp));

    host_grid2_complete g2(std::move(bin_data), std::move(axes));

    const auto& axis_r = g2.template get_axis<n_axis::label::e_r>();
    const auto& axis_phi = g2.template get_axis<n_axis::label::e_phi>();

    auto width_r = axis_r.bin_width();
    auto width_phi = axis_phi.bin_width();

    // pre-check
    for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
        for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
            auto bin = g2.bin({i_r, i_phi});
            auto invalid_bin =
                populator<host_grid2_complete::populator_impl>::init<point3>(first_tp);

            test_entry_collection(bin, invalid_bin.view());
        }
    }

    // run device side test, which populates the grid
    grid_complete_test(get_data(g2), axis_r.nbins(), axis_phi.nbins());

    // post-check
    for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
        for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
            const dindex gbin_idx = g2.serialize(
                detray::n_axis::multi_bin<2>{i_r, i_phi});
            const auto& bin = g2.bin(gbin_idx);

            // Other point with which the bin has been completed
            const detray::scalar
gbin_idx_f{static_cast<detray::scalar>(gbin_idx)}; const point3 tp{axis_r.min()
+ gbin_idx_f * width_r, axis_phi.min() + gbin_idx_f * width_phi, 0.5};

            // Go through all points and compare
            int pt_idx{0};
            for (const auto& pt : bin) {
                if (pt_idx == 0) {
                    EXPECT_EQ(pt, first_tp);
                } else {
                    EXPECT_EQ(pt, tp);
                }
                pt_idx++;
            }
        }
    }
}

// Test both the attach population and the const grid reading
TEST(grids_cuda, grid2_attach_populator) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // Build multi-axis
    using axes_t = polar<host_container_types>;

    typename axes_t::boundary_storage_type axis_data(&mng_mr);
    typename axes_t::edges_storage_type bin_edges(&mng_mr);

    axis_data.reserve(2);
    axis_data.insert(axis_data.begin(),
                      {dindex_range{0u, 2u}, dindex_range{2u, 65u}});
    bin_edges.reserve(4);
    bin_edges.insert(bin_edges.begin(), {0.f, 6.f, -constant<scalar>::pi,
constant<scalar>::pi});

    axes_t axes(std::move(axis_data), std::move(bin_edges));

    // build host grid
    const point3 first_tp{3.f, 3.f, 3.f};
    const point3 invalid_tp{0.f, 0.f, 0.f};

    host_grid2_attach::bin_storage_type bin_data(&mng_mr);
    bin_data.resize(2u * 65u,
                    populator<host_grid2_attach::populator_impl>::init<point3>(first_tp));

    host_grid2_attach g2(std::move(bin_data), std::move(axes));

    const auto& axis_r = g2.template get_axis<n_axis::label::e_r>();
    const auto& axis_phi = g2.template get_axis<n_axis::label::e_phi>();

    auto width_r = axis_r.bin_width();
    auto width_phi = axis_phi.bin_width();

    // pre-check
    for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
        for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
            auto bin = g2.bin({i_r, i_phi});
            auto invalid_bin =
                populator<host_grid2_complete::populator_impl>::init<point3>(first_tp);

            test_entry_collection(bin, invalid_bin.view());
        }
    }

    // run device side test, which populates the grid
    grid_attach_test(get_data(g2), axis_r.nbins(), axis_phi.nbins());

    // post-check
    for (unsigned int i_r = 0u; i_r < axis_r.nbins(); i_r++) {
        for (unsigned int i_phi = 0u; i_phi < axis_phi.nbins(); i_phi++) {
            const dindex gbin_idx = g2.serialize(
                detray::n_axis::multi_bin<2>{i_r, i_phi});
            const auto& bin = g2.bin(gbin_idx);

            // Other point with which the bin has been completed
            const detray::scalar
gbin_idx_f{static_cast<detray::scalar>(gbin_idx)}; const point3 tp{axis_r.min()
+ gbin_idx_f * width_r, axis_phi.min() + gbin_idx_f * width_phi, 0.5};

            // Go through all points and compare
            int pt_idx{0};
            for (const auto& pt : bin) {
                if (pt_idx == 0) {
                    EXPECT_EQ(pt, first_tp);
                } else if (pt_idx == 1) {
                    EXPECT_EQ(pt, tp);
                } else {
                    EXPECT_EQ(pt, invalid_tp);
                }
                pt_idx++;
            }
        }
    }

    // This seems to be undefined behaviour for now
    //const host_grid2_attach& g2_const = g2;
    // Read the grid as a const grid of non-const values
    //grid_attach_read_test(get_data(g2_const), axis_r.nbins(),
axis_phi.nbins());
}*/

TEST(grids_cuda, cylindrical3D_collection) {
    // Data-owning grid collection
    vecmem::cuda::managed_memory_resource mng_mr;

    using grid_collection_t = grid_collection<n_own_host_grid2_attach>;

    vecmem::vector<typename grid_collection_t::size_type> grid_offsets(&mng_mr);
    typename grid_collection_t::bin_storage_type bin_data(&mng_mr);
    typename grid_collection_t::axes_storage_type edge_ranges(&mng_mr);
    typename grid_collection_t::edges_storage_type bin_edges(&mng_mr);

    // Offsets for the grids into the bin storage
    grid_offsets.reserve(3u);
    grid_offsets.insert(grid_offsets.begin(), {0u, 48u, 72u});

    // Offsets into edges container and #bins for all axes
    edge_ranges.reserve(9u);
    edge_ranges.insert(
        edge_ranges.begin(),
        {dindex_range{0u, 2u}, dindex_range{2u, 4u}, dindex_range{4u, 6u},
         dindex_range{6u, 1u}, dindex_range{8u, 3u}, dindex_range{10u, 8u},
         dindex_range{12u, 5u}, dindex_range{14u, 5u}, dindex_range{16u, 5u}});

    // Bin edges for all axes (two boundaries for regular binned axes)
    bin_edges.reserve(18u);
    bin_edges.insert(bin_edges.begin(),
                     {-10.f, 10.f, -20.f, 20.f, 0.f, 120.f, -5.f, 5.f, -15.f,
                      15.f, 0.f, 50.f, -15.f, 15.f, -35.f, 35.f, 0.f, 550.f});

    // Bin test entries
    bin_data.resize(197u);
    std::generate_n(
        bin_data.begin(), 197u,
        bin_content_sequence<populator<n_own_host_grid2_attach::populator_impl>,
                             dindex>());

    vecmem::vector<unsigned int> n_bins(9u, &mng_mr);
    vecmem::vector<std::array<dindex, 3>> result_bins(bin_data.size(), &mng_mr);

    grid_collection_t grid_coll(std::move(grid_offsets), std::move(bin_data),
                                std::move(edge_ranges), std::move(bin_edges));

    // Call test function
    const auto& axis_r = grid_coll[2].template get_axis<n_axis::label::e_r>();
    const auto& axis_phi =
        grid_coll[2].template get_axis<n_axis::label::e_phi>();
    const auto& axis_z = grid_coll[2].template get_axis<n_axis::label::e_z>();

    grid_collection_test(get_data(grid_coll), vecmem::get_data(n_bins),
                         vecmem::get_data(result_bins), grid_coll.size(),
                         axis_r.nbins(), axis_phi.nbins(), axis_z.nbins());

    // Compare results
    EXPECT_EQ(4u, n_bins[0]);
    EXPECT_EQ(4u, n_bins[1]);
    EXPECT_EQ(8u, n_bins[2]);
    EXPECT_EQ(3u, n_bins[3]);
    EXPECT_EQ(3u, n_bins[4]);
    EXPECT_EQ(10u, n_bins[5]);
    EXPECT_EQ(7u, n_bins[6]);
    EXPECT_EQ(5u, n_bins[7]);
    EXPECT_EQ(7u, n_bins[8]);

    for (unsigned int i{0u}; i < bin_data.size(); ++i) {
        EXPECT_EQ(bin_data[i].content(), result_bins[i]);
    }
}
