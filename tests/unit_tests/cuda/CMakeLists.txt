add_detray_test(grids_grid2_cuda
  "grids_grid2_cuda.cpp"
  "grids_grid2_cuda_kernel.hpp"
  "grids_grid2_cuda_kernel.cu"
  LINK_LIBRARIES
  vecmem::cuda)

add_detray_test(transform_store
  "transform_store_cuda.cpp"
  "transform_store_cuda_kernel.hpp"
  "transform_store_cuda_kernel.cu"
  LINK_LIBRARIES
  vecmem::cuda detray::array detray::eigen detray::smatrix)
