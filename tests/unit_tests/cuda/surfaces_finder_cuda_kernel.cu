/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/cuda_definitions.hpp"
#include "surfaces_finder_cuda_kernel.hpp"

namespace detray {

__global__ void surfaces_finder_test_kernel(
    surfaces_finder_view<surfaces_finder_host_t> finder_data,
    vecmem::data::vector_view<dindex> outputs_data) {

    surfaces_finder_device_t finder_device(finder_data);
    vecmem::device_vector<dindex> outputs_device(outputs_data);

    // fill the output vector with grid elements
    int counts = 0;
    for (unsigned int i_g = 0; i_g < finder_device.size(); i_g++) {
        auto& g2 = finder_device[i_g];
        auto& xaxis = g2.axis_p0();
        auto& yaxis = g2.axis_p1();

        for (unsigned int i_y = 0; i_y < yaxis.bins(); i_y++) {
            for (unsigned int i_x = 0; i_x < xaxis.bins(); i_x++) {

                auto data = g2.bin(i_x, i_y);

                for (auto& d : data) {
                    outputs_device[counts] = d;
                    counts++;
                }
            }
        }
    }
}

void surfaces_finder_test(
    surfaces_finder_view<surfaces_finder_host_t> finder_data,
    vecmem::data::vector_view<dindex>& outputs_data) {

    constexpr int block_dim = 1;
    constexpr int thread_dim = 1;

    surfaces_finder_test_kernel<<<block_dim, thread_dim>>>(finder_data,
                                                           outputs_data);

    // cuda error check
    DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
    DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray