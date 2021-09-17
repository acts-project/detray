add_library(cuda_test SHARED
  grids_grid2_cuda_kernel.cu)

target_include_directories(cuda_test
  PUBLIC
  ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
  )

target_compile_options(cuda_test
  PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:
  #-arch=sm_52 
  #-gencode=arch=compute_75,code=sm_75 
  --expt-relaxed-constexpr>)
  #--use_fast_math>)
  #-fmad=false>)

target_link_libraries(cuda_test PUBLIC detray::core vecmem::cuda detray::tests_common)

add_detray_test(grids_grid2_cuda grids_grid2_cuda.cpp "detray::core;vecmem::cuda;cuda_test")
