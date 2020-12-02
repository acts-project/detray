add_library(detray_tests_io INTERFACE)

target_include_directories(detray_tests_io
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

install(
  DIRECTORY include/io
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tests)

add_library(detray::tests_io ALIAS detray_tests_io)