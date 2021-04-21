add_library(detray_apps_common INTERFACE)

target_include_directories(detray_apps_common
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

add_library(detray::apps ALIAS detray_apps_common)