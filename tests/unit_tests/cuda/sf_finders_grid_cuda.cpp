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
    dindex i = 0;
    for (const auto& elem : bin_content) {
        ASSERT_FLOAT_EQ(elem, expected[i++]);
    }
}

/// Test collection entries in a bin element by element
template <typename content_t, typename expected_content_t>
void test_entry_collection(const content_t& bin_content,
                           const expected_content_t& expected) {

    // Running indices for the points in the collection and the elements of a
    // single point3
    dindex i = 0, j = 0;
    for (const auto& entry : bin_content) {
        const auto& expt_entry = expected[i++];
        j = 0;
        for (const auto& elem : entry) {
            ASSERT_FLOAT_EQ(elem, expt_entry[j++]);
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
    axis_data.insert(axis_data.begin(), {dindex_range{0, 3}, dindex_range{2, 6},
                                         dindex_range{4, 10}});
    bin_edges.reserve(6);
    bin_edges.insert(bin_edges.begin(), {-1., 4., 0., 6., -5., 5.});

    axes_t axes(std::move(axis_data), std::move(bin_edges));

    // build host grid
    host_grid3_replace::bin_storage_type bin_data(&mng_mr);
    bin_data.resize(3 * 6 * 10,
                    host_grid3_replace::populator_type::init<point3>());

    host_grid3_replace g3(std::move(bin_data), std::move(axes));

    const auto& axis_x = g3.template get_axis<n_axis::label::e_x>();
    const auto& axis_y = g3.template get_axis<n_axis::label::e_y>();
    const auto& axis_z = g3.template get_axis<n_axis::label::e_z>();

    // pre-check
    for (std::size_t i_x = 0; i_x < axis_x.nbins(); i_x++) {
        for (std::size_t i_y = 0; i_y < axis_y.nbins(); i_y++) {
            for (std::size_t i_z = 0; i_z < axis_z.nbins(); i_z++) {
                const auto& bin = g3.at({i_x, i_y, i_z});
                auto invalid_bin =
                    host_grid3_replace::populator_type::init<point3>();
                test_content(bin[0], invalid_bin.content());
            }
        }
    }
    // run device side test, which populates the grid
    grid_replace_test(get_data(g3), axis_x.nbins(), axis_y.nbins(),
                      axis_z.nbins());
    // post-check
    for (std::size_t i_x = 0; i_x < axis_x.nbins(); i_x++) {
        for (std::size_t i_y = 0; i_y < axis_y.nbins(); i_y++) {
            for (std::size_t i_z = 0; i_z < axis_z.nbins(); i_z++) {
                dindex gbin_idx = g3.serializer()(
                    g3.axes(), detray::n_axis::multi_bin<3>{i_x, i_y, i_z});

                const auto& bin = g3.at(gbin_idx);

                const point3 tp{axis_x.min() + gbin_idx * axis_x.bin_width(),
                                axis_y.min() + gbin_idx * axis_y.bin_width(),
                                axis_z.min() + gbin_idx * axis_z.bin_width()};

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

    axis_data.reserve(2);
    axis_data.insert(axis_data.end(), {dindex_range{0, 4}, dindex_range{5, 2}});
    bin_edges.reserve(7);
    bin_edges.insert(bin_edges.end(), {1., 3., 9., 27., 81., -2., 4.});

    axes_t axes(std::move(axis_data), std::move(bin_edges));

    // build host grid
    host_grid2_replace_ci::bin_storage_type bin_data(&mng_mr);
    bin_data.resize(2 * 4,
                    host_grid2_replace_ci::populator_type::init<point3>());

    host_grid2_replace_ci g2(std::move(bin_data), std::move(axes));

    const auto& axis_r = g2.template get_axis<n_axis::label::e_r>();
    const auto& axis_phi = g2.template get_axis<n_axis::label::e_phi>();

    // pre-check
    for (std::size_t i_r = 0; i_r < axis_r.nbins(); i_r++) {
        for (std::size_t i_phi = 0; i_phi < axis_phi.nbins(); i_phi++) {
            const auto& bin = g2.at({i_r, i_phi});
            auto invalid_bin =
                host_grid2_replace_ci::populator_type::init<point3>();

            test_content(bin[0], invalid_bin.content());
        }
    }

    // run device side test, which populates the grid
    grid_replace_ci_test(get_data(g2), axis_r.nbins(), axis_phi.nbins());

    // post-check
    for (std::size_t i_r = 0; i_r < axis_r.nbins(); i_r++) {
        for (std::size_t i_phi = 0; i_phi < axis_phi.nbins(); i_phi++) {
            dindex gbin_idx = g2.serializer()(
                g2.axes(), detray::n_axis::multi_bin<2>{i_r, i_phi});
            const auto& bin = g2.at(gbin_idx);

            const point3 tp{axis_r.min() + gbin_idx * axis_r.bin_width(i_r),
                            axis_phi.min() + gbin_idx * axis_phi.bin_width(),
                            0.5};

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

    axis_data.reserve(2);
    axis_data.insert(axis_data.begin(),
                      {dindex_range{0, 3}, dindex_range{2, 7}});
    bin_edges.reserve(4);
    bin_edges.insert(bin_edges.begin(), {0., 3., -1, 6.});

    axes_t axes(std::move(axis_data), std::move(bin_edges));

    // build host grid
    const point3 first_tp{3., 3., 3.};

    host_grid2_complete::bin_storage_type bin_data(&mng_mr);
    bin_data.resize(
        3 * 7, host_grid2_complete::populator_type::init<point3>(first_tp));

    host_grid2_complete g2(std::move(bin_data), std::move(axes));

    const auto& axis_r = g2.template get_axis<n_axis::label::e_r>();
    const auto& axis_phi = g2.template get_axis<n_axis::label::e_phi>();

    auto width_r = axis_r.bin_width();
    auto width_phi = axis_phi.bin_width();

    // pre-check
    for (std::size_t i_r = 0; i_r < axis_r.nbins(); i_r++) {
        for (std::size_t i_phi = 0; i_phi < axis_phi.nbins(); i_phi++) {
            auto bin = g2.at({i_r, i_phi});
            auto invalid_bin =
                host_grid2_complete::populator_type::init<point3>(first_tp);

            test_entry_collection(bin, invalid_bin.view());
        }
    }

    // run device side test, which populates the grid
    grid_complete_test(get_data(g2), axis_r.nbins(), axis_phi.nbins());

    // post-check
    for (std::size_t i_r = 0; i_r < axis_r.nbins(); i_r++) {
        for (std::size_t i_phi = 0; i_phi < axis_phi.nbins(); i_phi++) {
            dindex gbin_idx = g2.serializer()(
                g2.axes(), detray::n_axis::multi_bin<2>{i_r, i_phi});
            const auto& bin = g2.at(gbin_idx);

            // Other point with which the bin has been completed
            const point3 tp{axis_r.min() + gbin_idx * width_r,
                            axis_phi.min() + gbin_idx * width_phi, 0.5};

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
                      {dindex_range{0, 2}, dindex_range{2, 65}});
    bin_edges.reserve(4);
    bin_edges.insert(bin_edges.begin(), {0., 6., -M_PI, M_PI});

    axes_t axes(std::move(axis_data), std::move(bin_edges));

    // build host grid
    const point3 first_tp{3., 3., 3.};
    const point3 invalid_tp{0., 0., 0.};

    host_grid2_attach::bin_storage_type bin_data(&mng_mr);
    bin_data.resize(2 * 65,
                    host_grid2_attach::populator_type::init<point3>(first_tp));

    host_grid2_attach g2(std::move(bin_data), std::move(axes));

    const auto& axis_r = g2.template get_axis<n_axis::label::e_r>();
    const auto& axis_phi = g2.template get_axis<n_axis::label::e_phi>();

    auto width_r = axis_r.bin_width();
    auto width_phi = axis_phi.bin_width();

    // pre-check
    for (std::size_t i_r = 0; i_r < axis_r.nbins(); i_r++) {
        for (std::size_t i_phi = 0; i_phi < axis_phi.nbins(); i_phi++) {
            auto bin = g2.at({i_r, i_phi});
            auto invalid_bin =
                host_grid2_complete::populator_type::init<point3>(first_tp);

            test_entry_collection(bin, invalid_bin.view());
        }
    }

    // run device side test, which populates the grid
    grid_attach_test(get_data(g2), axis_r.nbins(), axis_phi.nbins());

    // post-check
    for (std::size_t i_r = 0; i_r < axis_r.nbins(); i_r++) {
        for (std::size_t i_phi = 0; i_phi < axis_phi.nbins(); i_phi++) {
            dindex gbin_idx = g2.serializer()(
                g2.axes(), detray::n_axis::multi_bin<2>{i_r, i_phi});
            const auto& bin = g2.at(gbin_idx);

            // Other point with which the bin has been completed
            const point3 tp{axis_r.min() + gbin_idx * width_r,
                            axis_phi.min() + gbin_idx * width_phi, 0.5};

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
    grid_offsets.reserve(3);
    grid_offsets.insert(grid_offsets.begin(), {0UL, 48UL, 72UL});

    // Offsets into edges container and #bins for all axes
    edge_ranges.reserve(9);
    edge_ranges.insert(
        edge_ranges.begin(),
        {dindex_range{0u, 2u}, dindex_range{2u, 4u}, dindex_range{4u, 6u},
         dindex_range{6u, 1u}, dindex_range{8u, 3u}, dindex_range{10u, 8u},
         dindex_range{12u, 5u}, dindex_range{14u, 5u}, dindex_range{16u, 5u}});

    // Bin edges for all axes (two boundaries for regular binned axes)
    bin_edges.reserve(18);
    bin_edges.insert(bin_edges.begin(),
                     {-10, 10., -20., 20., 0., 120., -5., 5., -15., 15., 0.,
                      50., -15, 15., -35., 35., 0., 550.});

    // Bin test entries
    bin_data.resize(197UL);
    std::generate_n(
        bin_data.begin(), 197UL,
        bin_content_sequence<n_own_host_grid2_attach::populator_type,
                             dindex>());

    vecmem::vector<std::size_t> n_bins(9, &mng_mr);
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
    EXPECT_EQ(2UL, n_bins[0]);
    EXPECT_EQ(4UL, n_bins[1]);
    EXPECT_EQ(6UL, n_bins[2]);
    EXPECT_EQ(1UL, n_bins[3]);
    EXPECT_EQ(3UL, n_bins[4]);
    EXPECT_EQ(8UL, n_bins[5]);
    EXPECT_EQ(5UL, n_bins[6]);
    EXPECT_EQ(5UL, n_bins[7]);
    EXPECT_EQ(5UL, n_bins[8]);

    for (std::size_t i{0}; i < bin_data.size(); ++i) {
        EXPECT_EQ(bin_data[i].content(), result_bins[i]);
    }
}