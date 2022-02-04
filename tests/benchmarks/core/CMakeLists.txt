# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2022 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

detray_add_executable( benchmark_grids "benchmark_grids.cpp"
   LINK_LIBRARIES benchmark::benchmark vecmem::core detray_tests_common
                  detray::core )
