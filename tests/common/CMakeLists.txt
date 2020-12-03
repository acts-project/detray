add_library(detray_tests_common INTERFACE)

target_include_directories(detray_tests_common
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if(DETRAY_CUSTOM_SCALARTYPE)
  target_compile_definitions(
    detray_eigen
    PUBLIC -DDETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE})
endif()

install(
  DIRECTORY include/common
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tests)

add_library(detray::tests_common ALIAS detray_tests_common)