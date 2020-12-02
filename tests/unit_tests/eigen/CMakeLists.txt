enable_testing()

add_detray_test(eigen_plugin
                eigen_plugin.cpp
                detray::eigen)

add_detray_test(eigen_core
                eigen_core.cpp
                detray::eigen)

add_detray_test(eigen_planar_intersection
                eigen_planar_intersection.cpp
                detray::eigen)

add_detray_test(eigen_read_detector
                eigen_read_detector.cpp
                detray::eigen)
