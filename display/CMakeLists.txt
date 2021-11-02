message(STATUS "Building 'detray::display' component")

add_library(detray_display INTERFACE)

target_include_directories(detray_display
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(detray_display INTERFACE detray::core Matplot++::matplot)

install(
  DIRECTORY include/style
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  DIRECTORY include/view
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(detray::display ALIAS detray_display)