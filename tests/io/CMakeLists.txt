add_library(detray_tests_io INTERFACE)

target_include_directories(detray_tests_io
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(DETRAY_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    detray_tests_io
    INTERFACE -DDETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE})
endif()

add_library(detray::tests_io ALIAS detray_tests_io)

if (DETRAY_IO_JSON)
    add_executable(json_write_detector json/write_detector.cpp)
    target_link_libraries(json_write_detector detray::core detray::io detray::tests_common detray::tests_io detray::array)
    add_executable(json_write_grid json/write_grid.cpp)
    target_link_libraries(json_write_grid detray::core detray::io detray::tests_common detray::tests_io detray::array)
endif()