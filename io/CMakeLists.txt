add_library(detray_io INTERFACE)

target_include_directories(detray_io
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(detray_io INTERFACE nlohmann_json)

install(TARGETS detray_io
        PUBLIC_HEADER
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/core
)

install(
  DIRECTORY include/json
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  add_library(detray::io ALIAS detray_io)
