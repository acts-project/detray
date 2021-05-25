enable_testing()

set(PLUGIN_EXTRA_LIBRARIES ROOT::MathCore)

foreach(etest ${all_unit_tests})
    add_detray_test(smatrix_${etest}
                    smatrix_${etest}.cpp 
                    detray::smatrix)
endforeach(etest)

if (DETRAY_IO_CSV)
    set(UNIT_TEST_EXTRA_LIBRARIES detray::io)
    add_detray_test(smatrix_read_detector
                    smatrix_read_detector.cpp
                    detray::smatrix)
endif()

if (DETRAY_DISPLAY)
    set(UNIT_TEST_EXTRA_LIBRARIES detray::display)
    add_detray_test(smatrix_display_masks
                    smatrix_display_masks.cpp
                    detray::smatrix)
endif()
