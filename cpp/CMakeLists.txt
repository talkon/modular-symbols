set(SOURCES manin_symbol.cpp manin_symbol.h cache_decorator.h)

add_executable(main main.cpp ${SOURCES})
add_library(modular_symbols STATIC ${SOURCES})