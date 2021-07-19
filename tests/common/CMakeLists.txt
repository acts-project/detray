add_library(detray_tests_common INTERFACE)

target_include_directories(detray_tests_common
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(DETRAY_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    detray_tests_common
    INTERFACE -DDETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE})
endif()

if(DETRAY_BENCHMARKS AND DETRAY_BENCHMARKS_REP)
  target_compile_definitions(
    detray_tests_common
    INTERFACE -DDETRAY_BENCHMARKS_REP=${DETRAY_BENCHMARKS_REP})
endif()

if(DETRAY_BENCHMARKS_MULTITHREAD)
  target_compile_definitions(
    detray_tests_common
    INTERFACE -DDETRAY_BENCHMARKS_MULTITHREAD=ON)
endif()

add_library(detray::tests_common ALIAS detray_tests_common)