enable_testing()

add_detray_test(smatrix_masks_container
                smatrix_masks_container.cpp detray::smatrix)

add_detray_test(smatrix_annulus2
                smatrix_annulus2.cpp detray::smatrix)

add_detray_test(smatrix_cylinder3
                smatrix_cylinder3.cpp detray::smatrix)

add_detray_test(smatrix_rectangle2
                smatrix_rectangle2.cpp detray::smatrix)

add_detray_test(smatrix_ring2
                smatrix_ring2.cpp detray::smatrix)

add_detray_test(smatrix_single3
                smatrix_single3.cpp detray::smatrix)

add_detray_test(smatrix_trapezoid2
                smatrix_trapezoid2.cpp detray::smatrix)

add_detray_test(smatrix_unmasked
                smatrix_unmasked.cpp detray::smatrix)

add_detray_test(smatrix_plugin
                smatrix_plugin.cpp
                detray::smatrix)

add_detray_test(smatrix_core
                smatrix_core.cpp
                detray::smatrix)

add_detray_test(smatrix_cylinder_intersection
                smatrix_cylinder_intersection.cpp
                detray::smatrix)

add_detray_test(smatrix_planar_intersection
                smatrix_planar_intersection.cpp
                detray::smatrix)

add_detray_test(smatrix_read_detector
                smatrix_read_detector.cpp
                detray::smatrix)
                
target_link_libraries(smatrix_masks_container ROOT::MathCore)
target_link_libraries(smatrix_annulus2 ROOT::MathCore)
target_link_libraries(smatrix_cylinder3 ROOT::MathCore)
target_link_libraries(smatrix_rectangle2 ROOT::MathCore)
target_link_libraries(smatrix_ring2 ROOT::MathCore)
target_link_libraries(smatrix_single3 ROOT::MathCore)
target_link_libraries(smatrix_trapezoid2 ROOT::MathCore)
target_link_libraries(smatrix_unmasked ROOT::MathCore)
target_link_libraries(smatrix_plugin ROOT::MathCore)
target_link_libraries(smatrix_core ROOT::MathCore)
target_link_libraries(smatrix_cylinder_intersection ROOT::MathCore)
target_link_libraries(smatrix_planar_intersection ROOT::MathCore)
target_link_libraries(smatrix_read_detector ROOT::MathCore)

