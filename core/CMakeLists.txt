message(STATUS "Building 'detray::core' component")

add_library(detray_core INTERFACE)

target_include_directories(detray_core
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
)

target_link_libraries(detray_core INTERFACE vecmem::core)

install(
  DIRECTORY include/core
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  DIRECTORY include/masks
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  DIRECTORY include/utils
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  DIRECTORY include/tools
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(detray::core ALIAS detray_core)

