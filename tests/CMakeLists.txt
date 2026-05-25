find_package(Catch2 REQUIRED)

add_executable(test_time_scales
    test_time_scales.cpp
)

target_link_libraries(test_time_scales PRIVATE
    times
    Catch2::Catch2WithMain
)

include(CTest)
include(Catch)
catch_discover_tests(test_time_scales)