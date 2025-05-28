/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detector_construction.hpp"

#include "detray/core/detail/alignment.hpp"
#include "detray/test/common/build_toy_detector.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// System include(s)
#include <iostream>

/// Prepare the data and move it to device
int main() {

    using algebra_t = detray::tutorial::algebra_t;
    using scalar = detray::tutorial::scalar;

    // memory resource(s)
    vecmem::host_memory_resource host_mr;
    vecmem::cuda::managed_memory_resource mng_mr;
    vecmem::cuda::device_memory_resource dev_mr;

    // Helper object for performing memory copies to CUDA devices
    vecmem::cuda::copy cuda_cpy;

    std::cout
        << "Detector Device Construction Tutorial\n====================\n\n";

    //
    // Managed Memory
    //

    // create toy geometry with vecmem managed memory resouce
    auto [det_mng, names_mng] = detray::build_toy_detector<algebra_t>(mng_mr);

    // Get the view onto the detector data directly
    auto det_mng_data = detray::get_data(det_mng);

    // Pass the view and call the kernel
    std::cout << "Using CUDA unified memory:" << std::endl;
    detray::tutorial::print(det_mng_data);

    //
    // Default copy to device
    //

    // create toy geometry in host memory
    auto [det_host, names_host] =
        detray::build_toy_detector<algebra_t>(host_mr);

    // Copy the detector data to device (synchronous copy, fixed size buffers)
    auto det_fixed_buff = detray::get_buffer(det_host, dev_mr, cuda_cpy);

    // Get the detector view from the buffer and call the kernel
    std::cout << "\nSynchronous copy, fixed size buffers:" << std::endl;
    detray::tutorial::print(detray::get_data(det_fixed_buff));

    // Copy the data to device in resizable buffers (synchronous copy)
    auto det_resz_buff =
        detray::get_buffer(det_host, dev_mr, cuda_cpy, detray::copy::sync,
                           vecmem::data::buffer_type::resizable);

    std::cout << "\nSynchronous copy, resizable buffers:" << std::endl;
    detray::tutorial::print(detray::get_data(det_resz_buff));

    //
    // Custom copy to device
    //

    // Get each buffer individually
    auto vol_buff = detray::get_buffer(det_host.volumes(), dev_mr, cuda_cpy,
                                       detray::copy::sync,
                                       vecmem::data::buffer_type::fixed_size);
    auto sf_buff = detray::get_buffer(det_host.surfaces(), dev_mr, cuda_cpy,
                                      detray::copy::sync,
                                      vecmem::data::buffer_type::fixed_size);
    auto trf_buff = detray::get_buffer(det_host.transform_store(), dev_mr,
                                       cuda_cpy, detray::copy::sync,
                                       vecmem::data::buffer_type::fixed_size);
    auto msk_buff = detray::get_buffer(det_host.mask_store(), dev_mr, cuda_cpy,
                                       detray::copy::sync,
                                       vecmem::data::buffer_type::fixed_size);
    auto mat_buff = detray::get_buffer(det_host.material_store(), dev_mr,
                                       cuda_cpy, detray::copy::sync,
                                       vecmem::data::buffer_type::fixed_size);
    auto acc_buff = detray::get_buffer(det_host.accelerator_store(), dev_mr,
                                       cuda_cpy, detray::copy::sync,
                                       vecmem::data::buffer_type::fixed_size);
    auto vgrid_buff = detray::get_buffer(det_host.volume_search_grid(), dev_mr,
                                         cuda_cpy, detray::copy::sync,
                                         vecmem::data::buffer_type::fixed_size);

    // Assemble the detector buffer
    using host_detector_type = decltype(det_host);
    auto det_custom_buff = typename host_detector_type::buffer_type(
        std::move(vol_buff), std::move(sf_buff), std::move(trf_buff),
        std::move(msk_buff), std::move(mat_buff), std::move(acc_buff),
        std::move(vgrid_buff));

    std::cout << "\nCustom buffer setup:" << std::endl;
    detray::tutorial::print(detray::get_data(det_custom_buff));

    // Construct an "aligned" transform store
    using host_transform_type =
        host_detector_type::transform_container::value_type;

    typename host_detector_type::transform_container host_aligned_transforms;
    detray::tutorial::point3 shift{.1f * detray::unit<scalar>::mm,
                                   .2f * detray::unit<scalar>::mm,
                                   .3f * detray::unit<scalar>::mm};

    for (const auto& tf : det_host.transform_store()) {
        detray::tutorial::point3 shifted{tf.translation()[0] + shift[0],
                                         tf.translation()[1] + shift[1],
                                         tf.translation()[2] + shift[2]};
        host_aligned_transforms.push_back(
            host_transform_type{shifted, tf.x(), tf.y(), tf.z()});
    }

    auto trf_buff_shifted = detray::get_buffer(
        host_aligned_transforms, dev_mr, cuda_cpy, detray::copy::sync,
        vecmem::data::buffer_type::fixed_size);

    auto detector_view =
        detray::detail::misaligned_detector_view<host_detector_type>(
            det_custom_buff, trf_buff_shifted);

    std::cout << "\nCustom buffer setup (shifted):" << std::endl;
    detray::tutorial::print(detector_view);
}
