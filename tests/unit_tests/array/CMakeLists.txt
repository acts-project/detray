enable_testing()

foreach(etest ${all_unit_tests})
    add_detray_test(array_${etest}
                    array_${etest}.cpp 
                    algebra::array)
endforeach(etest)

if (DETRAY_IO_CSV)
    set(UNIT_TEST_EXTRA_LIBRARIES detray::io)
    add_detray_test(array_read_detector
                    array_read_detector.cpp
                    algebra::array)
endif()

if (DETRAY_DISPLAY)
    set(UNIT_TEST_EXTRA_LIBRARIES detray::display)
    add_detray_test(array_display_masks
                    array_display_masks.cpp
                    algebra::array)
endif()
