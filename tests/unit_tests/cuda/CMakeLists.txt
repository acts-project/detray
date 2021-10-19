add_detray_test(grids_grid2_cuda
  "grids_grid2_cuda.cpp"
  "grids_grid2_cuda_kernel.hpp"
  "grids_grid2_cuda_kernel.cu"
  LINK_LIBRARIES
  vecmem::cuda)

list(APPEND algebras "array;eigen;smatrix")

foreach(algebra ${algebras})
  add_detray_test(${algebra}_transform_store_cuda
    "transform_store_cuda.cpp"
    "transform_store_cuda_kernel.hpp"
    "transform_store_cuda_kernel.cu"
    LINK_LIBRARIES
    vecmem::cuda detray::${algebra})
    target_compile_definitions(detray_test_${algebra}_transform_store_cuda PRIVATE ${algebra}=${algebra})
endforeach()

