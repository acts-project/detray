cmake_minimum_required(VERSION 3.11)

project(detray VERSION 0.1 LANGUAGES CXX)

option(CMAKE_WARN_DEPRECATED "Not warning for deprecation in cmake" Off)
enable_testing()

include(CMakeDependentOption)

option(DETRAY_ARRAY_PLUGIN "Build std::array math plugin" On)
option(DETRAY_EIGEN_PLUGIN "Build Eigen math plugin" On)
option(DETRAY_SMATRIX_PLUGIN "Build ROOT/SMatrix math plugin" Off)
option(DETRAY_VC_ALGEBRA_PLUGIN "Build vc based simd algebra math plugin" Off)
option(DETRAY_IO_CSV "Build CSV IO module" On)
option(DETRAY_IO_JSON "Build JSON IO module" On)
option(DETRAY_DISPLAY "Build matplot++ based display module" Off)
option(DETRAY_USE_VECMEM "Build with vecmem" On)

if(DETRAY_CUSTOM_SCALARTYPE)
    message(STATUS "Building with custom scalar type: " ${DETRAY_CUSTOM_SCALARTYPE})
    set(ALGEBRA_PLUGIN_CUSTOM_SCALARTYPE ${DETRAY_CUSTOM_SCALARTYPE})
endif()

if (DETRAY_USE_VECMEM)
  set(ALGEBRA_PLUGIN_USE_VECMEM ON)
endif()

if (DETRAY_ARRAY_PLUGIN)
    set(ALGEBRA_PLUGIN_INCLUDE_ARRAY ON)
endif()

if (DETRAY_EIGEN_PLUGIN)
    set(ALGEBRA_PLUGIN_INCLUDE_EIGEN ON)
    find_package(Eigen3 REQUIRED)
    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})
endif()

if (DETRAY_SMATRIX_PLUGIN)
    set(ALGEBRA_PLUGIN_INCLUDE_SMATRIX ON)
    find_package(ROOT REQUIRED COMPONENTS MathCore)
    include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})
endif()

if (DETRAY_VC_ALGEBRA_PLUGIN)
    set(ALGEBRA_PLUGIN_INCLUDE_VC ON)
endif()

option(DETRAY_UNIT_TESTS "Enable unit tests" On)
option(DETRAY_BENCHMARKS "Enable benchmark tests" On)
cmake_dependent_option(DETRAY_BENCHMARKS_MULTITHREAD "Enable multithreaded benchmarks" Off "DETRAY_BENCHMARKS" Off)
if(DETRAY_BENCHMARKS_MULTITHREAD)
    message(STATUS "Using multithreaded benchmarks")
endif()
if(DETRAY_BENCHMARKS AND DETRAY_BENCHMARKS_REP)
    message(STATUS "Repetitions for benchmarks: " ${DETRAY_BENCHMARKS_REP})
endif()

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")

add_subdirectory(apps)
add_subdirectory(core)
add_subdirectory(extern)
add_subdirectory(plugins)
add_subdirectory(io)
add_subdirectory(tests)
if (DETRAY_DISPLAY)
    add_subdirectory(display)
endif()
