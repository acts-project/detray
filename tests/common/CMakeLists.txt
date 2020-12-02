add_library(detray_tests_common INTERFACE)

target_include_directories(detray_tests_common
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

install(
  DIRECTORY include/common
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tests)

add_library(detray::tests_common ALIAS detray_tests_common)