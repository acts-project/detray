add_detray_test(grids_grid2_cuda
  "grids_grid2_cuda.cpp"
  "grids_grid2_cuda_kernel.hpp"
  "grids_grid2_cuda_kernel.cu"
  LINK_LIBRARIES
  vecmem::cuda)
