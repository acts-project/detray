add_library(detray_core INTERFACE)

target_include_directories(detray_core
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS detray_core
        PUBLIC_HEADER
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/core
)

install(
  DIRECTORY include/masks
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  DIRECTORY include/core
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(detray::core ALIAS detray_core)

