find_package(Catch2 REQUIRED)

add_executable(test_times
    test_times.cpp
)

target_link_libraries(test_times PRIVATE
    times
    Catch2::Catch2WithMain
)

include(CTest)
include(Catch)
catch_discover_tests(test_times)