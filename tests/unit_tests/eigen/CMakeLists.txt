enable_testing()

foreach(etest ${all_unit_tests})
    add_detray_test(eigen_${etest}
                eigen_${etest}.cpp detray::eigen)
endforeach(etest)

if (DETRAY_DISPLAY)
    set(UNIT_TEST_EXTRA_LIBRARIES detray::display)
    add_detray_test(eigen_display_masks
                    eigen_display_masks.cpp
                    detray::eigen)
endif()
