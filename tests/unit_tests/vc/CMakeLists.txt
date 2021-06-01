enable_testing()

foreach(etest ${all_unit_tests})
    add_detray_test(vc_array_${etest}
                    vc_array_${etest}.cpp 
                    detray::vc_array)
endforeach(etest)

if (DETRAY_IO_CSV)
    set(UNIT_TEST_EXTRA_LIBRARIES detray::io)
    add_detray_test(vc_array_read_detector
                    vc_array_read_detector.cpp
                    detray::vc_array)
endif()

if (DETRAY_DISPLAY)
    set(UNIT_TEST_EXTRA_LIBRARIES detray::display)
    add_detray_test(vc_array_display_masks
                    vc_array_display_masks.cpp
                    detray::vc_array)
endif()
