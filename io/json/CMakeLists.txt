add_library(detray_io_json INTERFACE)

target_include_directories(detray_io_json
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(detray_io_json INTERFACE nlohmann_json)

install(
  DIRECTORY include/io
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  add_library(detray::io_json ALIAS detray_io_json)
