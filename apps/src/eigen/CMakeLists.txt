if (DETRAY_DISPLAY)
    add_executable(eigen_display_volumes_rz eigen_display_volumes_rz.cpp)
    target_link_libraries(eigen_display_volumes_rz detray::core detray::apps detray::io detray::display detray::eigen)

    add_executable(eigen_display_layers eigen_display_layers.cpp)
    target_link_libraries(eigen_display_layers detray::core detray::apps detray::io detray::display detray::eigen)
endif()