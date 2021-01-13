enable_testing()

add_detray_test(eigen_masks_container
                eigen_masks_container.cpp detray::eigen)

add_detray_test(eigen_annulus
                eigen_annulus.cpp detray::eigen)

add_detray_test(eigen_cylindrical_detector
                eigen_cylindrical_detector.cpp detray::eigen)

add_detray_test(eigen_cylinder3
                eigen_cylinder3.cpp detray::eigen)

add_detray_test(eigen_rectangle2
                eigen_rectangle2.cpp detray::eigen)

add_detray_test(eigen_ring2
                eigen_ring2.cpp detray::eigen)

add_detray_test(eigen_single3
                eigen_single3.cpp detray::eigen)

add_detray_test(eigen_trapezoid2
                eigen_trapezoid2.cpp detray::eigen)

add_detray_test(eigen_unmasked
                eigen_unmasked.cpp detray::eigen)

add_detray_test(eigen_plugin
                eigen_plugin.cpp
                detray::eigen)

add_detray_test(eigen_core
                eigen_core.cpp
                detray::eigen)

add_detray_test(eigen_cylinder_intersection
                eigen_cylinder_intersection.cpp
                detray::eigen)

add_detray_test(eigen_planar_intersection
                eigen_planar_intersection.cpp
                detray::eigen)

add_detray_test(eigen_read_detector
                eigen_read_detector.cpp
                detray::eigen)
