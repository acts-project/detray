add_library(detray_io_csv INTERFACE)

target_include_directories(detray_io_csv
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(detray_io_csv INTERFACE nlohmann_csv)

install(
  DIRECTORY include/io
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  add_library(detray::io_csv ALIAS detray_io_csv)
