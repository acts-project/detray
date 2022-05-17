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
        detray_add_test( ${algebra}_vecpar_cpu_host_mem
                "host_memory/rk_stepper_vecpar_hm.cpp"
                "algorithm/rk_stepper_vecpar.hpp"
                LINK_LIBRARIES GTest::gtest_main vecmem::core detray_tests_common
                detray::${algebra} vecpar::all OpenMP::OpenMP_CXX)

        target_compile_definitions( detray_test_${algebra}_vecpar_cpu_host_mem
                PRIVATE ${algebra}=${algebra} )

        set_target_properties( detray_test_${algebra}_vecpar_cpu_host_mem
                PROPERTIES CXX_STANDARD 20 )
    endforeach()
endif()

# if CUDA libraries available
if( DETRAY_BUILD_CUDA )
    find_package(CUDAToolkit)
    ### CPU tests
    if ((TARGET vecpar::all) AND (TARGET vecpar::omp))
        foreach(algebra ${algebras})
            # Unit tests for the selected algebra.
            detray_add_test( ${algebra}_vecpar_cpu
                    "managed_memory/rk_stepper_vecpar_mm.cpp"
                    "host_memory/rk_stepper_vecpar_hm.cpp"
                    "algorithm/rk_stepper_vecpar.hpp"
                    LINK_LIBRARIES GTest::gtest_main vecmem::cuda detray_tests_common
                    detray::${algebra} vecpar::all CUDA::cudart OpenMP::OpenMP_CXX)

            target_compile_definitions( detray_test_${algebra}_vecpar_cpu
                    PRIVATE ${algebra}=${algebra} )

            set_target_properties( detray_test_${algebra}_vecpar_cpu
                    PROPERTIES CXX_STANDARD 20 )

        endforeach()
    endif()

    ### GPU tests can be built only by Clang compiler +nvptx support
    if ((TARGET vecpar::all) AND (TARGET vecpar::cuda))
        foreach(algebra ${algebras})

            # Unit tests for the selected algebra.
            detray_add_test(${algebra}_vecpar_gpu
                    "managed_memory/rk_stepper_vecpar_mm.cpp"
             #       "host_memory/rk_stepper_vecpar_hm.cpp"
                    "algorithm/rk_stepper_vecpar.hpp"
                    LINK_LIBRARIES GTest::gtest_main vecmem::cuda detray_tests_common
                    detray::${algebra} vecpar::all CUDA::cudart)

            target_compile_definitions( detray_test_${algebra}_vecpar_gpu
                    PRIVATE ${algebra}=${algebra} )

            set_target_properties( detray_test_${algebra}_vecpar_gpu
                    PROPERTIES CXX_STANDARD 20 )

            target_compile_options( detray_test_${algebra}_vecpar_gpu PUBLIC
                    $<$<COMPILE_LANGUAGE:CXX>:-x cuda --offload-arch=sm_${CMAKE_CUDA_ARCHITECTURES}>)
        endforeach()
    endif()

endif() # if CUDA support enabled

