# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Enable CUDA as a language.
enable_language(CUDA)

add_detray_test(grids_grid2_cuda
  "grids_grid2_cuda.cpp"
  "grids_grid2_cuda_kernel.hpp"
  "grids_grid2_cuda_kernel.cu"
  LINK_LIBRARIES
  vecmem::cuda)

# make unit tests for multiple algebras
# Currently vc and smatrix is not supported
list(APPEND algebras "array;eigen")

foreach(algebra ${algebras})

  # transform store test
  add_detray_test(${algebra}_transform_store_cuda
    "transform_store_cuda.cpp"
    "transform_store_cuda_kernel.hpp"
    "transform_store_cuda_kernel.cu"
    LINK_LIBRARIES
    vecmem::cuda detray::${algebra})
  target_compile_definitions(detray_test_${algebra}_transform_store_cuda PRIVATE ${algebra}=${algebra})

  # mask store test
  add_detray_test(${algebra}_mask_store_cuda
    "mask_store_cuda.cpp"
    "mask_store_cuda_kernel.hpp"
    "mask_store_cuda_kernel.cu"
    LINK_LIBRARIES
    vecmem::cuda detray::${algebra})
  target_compile_definitions(detray_test_${algebra}_mask_store_cuda PRIVATE ${algebra}=${algebra})
  
  # tuple test
  add_detray_test(${algebra}_tuple_test_cuda
    "tuple_test_cuda.cpp"
    "tuple_test_cuda_kernel.hpp"
    "tuple_test_cuda_kernel.cu"
    LINK_LIBRARIES
    vecmem::cuda detray::${algebra})
  target_compile_definitions(detray_test_${algebra}_tuple_test_cuda PRIVATE ${algebra}=${algebra})

endforeach()

