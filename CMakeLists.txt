cmake_minimum_required(VERSION 3.11)

# must be set before project(...) call; version module is needed before
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(detray VERSION 0.0 LANGUAGES CXX)

include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})

option(DETRAY_EIGEN_PLUGIN "Build Eigen math plugin" On)
option(DETRAY_SMATRIX_PLUGIN "Build ROOT/SMatrix math plugin" Off)

option(DETRAY_TESTS "Enable unit tests" On)

if(DETRAY_EIGEN_PLUGIN)
    find_package(Eigen3 REQUIRED)
    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})
endif()

if(DETRAY_SMATRIX_PLUGIN)
    find_package(ROOT REQUIRED COMPONENTS MathCore)
    include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})
endif()

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)

add_subdirectory(core)
add_subdirectory(plugins)

if(DETRAY_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
