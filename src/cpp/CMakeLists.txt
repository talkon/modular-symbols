find_library(FLINT_LIBRARY flint REQUIRED)
find_library(GMP_LIBRARY gmp REQUIRED)
find_library(MPFR_LIBRARY mpfr REQUIRED)

add_executable(main manin_symbol.cpp manin_symbol.h cache_decorator.h)

target_link_libraries(main ${FLINT_LIBRARY} ${GMP_LIBRARY} ${MPFR_LIBRARY})