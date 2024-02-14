set(SOURCES manin_symbol.cpp manin_basis.cpp)

add_executable(main main.cpp ${SOURCES})
add_library(modular_symbols STATIC ${SOURCES})