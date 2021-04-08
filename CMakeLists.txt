cmake_minimum_required(VERSION 3.11)

project(detray VERSION 0.1 LANGUAGES CXX)

option(CMAKE_WARN_DEPRECATED "Not warning for deprecation in cmake" Off)
enable_testing()

option(DETRAY_ARRAY_PLUGIN "Build std::array math plugin" On)
option(DETRAY_EIGEN_PLUGIN "Build Eigen math plugin" On)
option(DETRAY_SMATRIX_PLUGIN "Build ROOT/SMatrix math plugin" Off)
option(DETRAY_IO_CSV "Build CSV IO module" On)
option(DETRAY_IO_JSON "Build JSON IO module" On)
option(DETRAY_DISPLAY "Build matplot++ based display module" On)

if(DETRAY_CUSTOM_SCALARTYPE)
    message(STATUS "Building with custom scalar type: " ${ACTS_CUSTOM_SCALARTYPE})
endif()

if (DETRAY_EIGEN_PLUGIN)
    find_package(Eigen3 REQUIRED)
    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})
endif()

if (DETRAY_SMATRIX_PLUGIN)
    find_package(ROOT REQUIRED COMPONENTS MathCore)
    include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})
endif()

option(DETRAY_UNIT_TESTS "Enable unit tests" On)
option(DETRAY_BENCHMARKS "Enable benchmark tests" On)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")

add_subdirectory(core)
add_subdirectory(extern)
add_subdirectory(plugins)
add_subdirectory(io)
add_subdirectory(tests)
if (DETRAY_DISPLAY)
    add_subdirectory(display)
endif()
