include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(tests manin_basis_test.cpp boundary_map_test.cpp newspace_test.cpp)
target_link_libraries(tests GTest::gtest_main modular_symbols)

include(GoogleTest)
gtest_discover_tests(tests)