enable_testing()

foreach(etest ${all_unit_tests})
    add_detray_test(array_${etest}
                    array_${etest}.cpp
                    LINK_LIBRARIES
                    detray::array)
endforeach(etest)

if (DETRAY_IO_CSV)
    add_detray_test(array_read_detector
                    array_read_detector.cpp
                    LINK_LIBRARIES
                    detray::array)
endif()

if (DETRAY_DISPLAY)
    add_detray_test(array_display_masks
                    array_display_masks.cpp
                    LINK_LIBRARIES
                    detray::array detray::display)
endif()
