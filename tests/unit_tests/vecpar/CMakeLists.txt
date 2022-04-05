# Detray library, part of the ACTS project (R&D line)
#
# (c) 2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0


# make unit tests for multiple algebras
set( algebras "array" )
if( DETRAY_EIGEN_PLUGIN )
    list( APPEND algebras "eigen" )
endif()

### CPU tests
find_package(OpenMP QUIET)
if ((TARGET vecpar::all) AND (TARGET vecpar::omp))
    foreach(algebra ${algebras})
       # Unit tests for the selected algebra.
       detray_add_test( ${algebra}_vecpar_cpu
               "rk_stepper_vecpar.cpp" "rk_stepper_vecpar.hpp"
               LINK_LIBRARIES GTest::gtest_main vecmem::cuda detray_tests_common
               detray::${algebra} vecpar::all OpenMP::OpenMP_CXX)

       target_compile_definitions( detray_test_${algebra}_vecpar_cpu
               PRIVATE ${algebra}=${algebra} )

       target_compile_options( detray_test_${algebra}_vecpar_cpu PUBLIC
               $<$<COMPILE_LANGUAGE:CXX>:-std=c++20>)
    endforeach()
endif()


### GPU tests can be built only by Clang compiler +nvptx support
find_package(CUDAToolkit)
if ((TARGET vecpar::all) AND (TARGET vecpar::cuda))
    foreach(algebra ${algebras})

        # Unit tests for the selected algebra.
        detray_add_test(${algebra}_vecpar_gpu
                "rk_stepper_vecpar.cpp" "rk_stepper_vecpar.hpp"
                LINK_LIBRARIES GTest::gtest_main vecmem::cuda detray_tests_common
                detray::${algebra} vecpar::all CUDA::cudart)

        target_compile_definitions( detray_test_${algebra}_vecpar_gpu
                PRIVATE ${algebra}=${algebra} )

        target_compile_options( detray_test_${algebra}_vecpar_gpu PUBLIC
                $<$<COMPILE_LANGUAGE:CXX>:-std=c++20>
                $<$<COMPILE_LANGUAGE:CXX>:-x cuda --offload-arch=sm_${CMAKE_CUDA_ARCHITECTURES}>)
    endforeach()
endif()