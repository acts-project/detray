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
#include <iostream>
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

}  // anonymous namespace

TEST(grids_cuda, grid3_replace_populator) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // Build multi-axis
    cartesian_3D<host_container_types> axes(mng_mr);

    auto axis_data = axes.m_data.axes_data();
    auto bin_edges = axes.m_data.edges();

    axis_data->reserve(3);
    axis_data->insert(
        axis_data->begin(),
        {dindex_range{0, 3}, dindex_range{2, 6}, dindex_range{4, 10}});
    bin_edges->reserve(6);
    bin_edges->insert(bin_edges->begin(), {-1., 4., 0., 6., -5., 5.});

    // build host grid
    host_grid3_replace::bin_storage_type bin_data(&mng_mr);
    bin_data.resize(3 * 6 * 10,
                    host_grid3_replace::populator_type::init<point3>());

    host_grid3_replace g3(std::move(bin_data), std::move(axes));

    const auto& axis_x = g3.template get_axis<n_axis::label::e_x>();
    const auto& axis_y = g3.template get_axis<n_axis::label::e_y>();
    const auto& axis_z = g3.template get_axis<n_axis::label::e_z>();

    auto width_x = axis_x.m_binning.bin_width();
    auto width_y = axis_y.m_binning.bin_width();
    auto width_z = axis_z.m_binning.bin_width();

    // pre-check
    for (std::size_t i_x = 0; i_x < axis_x.nbins(); i_x++) {
        for (std::size_t i_y = 0; i_y < axis_y.nbins(); i_y++) {
            for (std::size_t i_z = 0; i_z < axis_z.nbins(); i_z++) {
                const auto& bin = g3.at({i_x, i_y, i_z});
                auto invalid_bin =
                    host_grid3_replace::populator_type::init<point3>();
                test_content(*bin, invalid_bin.content());
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

                const point3 tp{axis_x.min() + gbin_idx * width_x,
                                axis_y.min() + gbin_idx * width_y,
                                axis_z.min() + gbin_idx * width_z};

                test_content(*bin, tp);
            }
        }
    }
}

TEST(grids_cuda, grid2_replace_populator_ci) {

    vecmem::cuda::managed_memory_resource mng_mr;

    // Build multi-axis
    polar_ir<host_container_types> axes(mng_mr);

    auto axis_data = axes.m_data.axes_data();
    auto bin_edges = axes.m_data.edges();

    axis_data->reserve(2);
    axis_data->insert(axis_data->end(),
                      {dindex_range{0, 4}, dindex_range{5, 2}});
    bin_edges->reserve(7);
    bin_edges->insert(bin_edges->end(), {1., 3., 9., 27., 81., -2., 4.});

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

            test_content(*bin, invalid_bin.content());
        }
    }

    // run device side test, which populates the grid
    grid_replace_ci_test(get_data(g2), axis_r.nbins(), axis_phi.nbins());

    // post-check
    auto width_phi = axis_phi.m_binning.bin_width();

    for (std::size_t i_r = 0; i_r < axis_r.nbins(); i_r++) {
        for (std::size_t i_phi = 0; i_phi < axis_phi.nbins(); i_phi++) {
            dindex gbin_idx = g2.serializer()(
                g2.axes(), detray::n_axis::multi_bin<2>{i_r, i_phi});
            const auto& bin = g2.at(gbin_idx);

            auto width_r = axis_r.m_binning.bin_width(i_r);

            const point3 tp{axis_r.min() + gbin_idx * width_r,
                            axis_phi.min() + gbin_idx * width_phi, 0.5};

            test_content(*bin, tp);
        }
    }
}

/*TEST(grids_cuda, grid2_complete_populator) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // Build multi-axis
    polar<host_container_types> axes(mng_mr);

    auto axis_data = axes.m_data.axes_data();
    auto bin_edges = axes.m_data.edges();

    axis_data->reserve(2);
    axis_data->insert(axis_data->begin(),
                      {dindex_range{0, 3}, dindex_range{2, 7}});
    bin_edges->reserve(4);
    bin_edges->insert(bin_edges->begin(), {0., 3., -1, 6.});

    // build host grid
    const point3 first_tp{3., 3., 3.};

    host_grid2_complete::bin_storage_type bin_data(&mng_mr);
    bin_data.resize(
        3 * 7, host_grid2_complete::populator_type::init<point3>(first_tp));

    host_grid2_complete g2(std::move(bin_data), std::move(axes));

    const auto& axis_r = g2.template get_axis<n_axis::label::e_r>();
    const auto& axis_phi = g2.template get_axis<n_axis::label::e_phi>();

    auto width_r = axis_r.m_binning.bin_width();
    auto width_phi = axis_phi.m_binning.bin_width();

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
    polar<host_container_types> axes(mng_mr);

    auto axis_data = axes.m_data.axes_data();
    auto bin_edges = axes.m_data.edges();

    axis_data->reserve(2);
    axis_data->insert(axis_data->begin(),
                      {dindex_range{0, 2}, dindex_range{2, 65}});
    bin_edges->reserve(4);
    bin_edges->insert(bin_edges->begin(), {0., 6., -M_PI, M_PI});

    // build host grid
    const point3 first_tp{3., 3., 3.};
    const point3 invalid_tp{0., 0., 0.};

    host_grid2_attach::bin_storage_type bin_data(&mng_mr);
    bin_data.resize(2 * 65,
                    host_grid2_attach::populator_type::init<point3>(first_tp));

    host_grid2_attach g2(std::move(bin_data), std::move(axes));

    const auto& axis_r = g2.template get_axis<n_axis::label::e_r>();
    const auto& axis_phi = g2.template get_axis<n_axis::label::e_phi>();

    auto width_r = axis_r.m_binning.bin_width();
    auto width_phi = axis_phi.m_binning.bin_width();

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
