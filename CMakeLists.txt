cmake_minimum_required(VERSION 3.14)
project(modular_symbols)

set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
set(CMAKE_CXX_FLAGS "-march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -ggdb")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g -fsanitize=address -fno-omit-frame-pointer")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_DEBUG FALSE)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# libraries required by all targets
find_library(FLINT_LIBRARY flint REQUIRED)
find_library(GMP_LIBRARY gmp REQUIRED)
find_library(MPFR_LIBRARY mpfr REQUIRED)

link_libraries(${FLINT_LIBRARY} ${GMP_LIBRARY} ${MPFR_LIBRARY})

enable_testing()

add_subdirectory(cpp)
add_subdirectory(test)