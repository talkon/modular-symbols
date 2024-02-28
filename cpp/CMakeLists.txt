set(SOURCES manin_symbol.cpp manin_basis.cpp boundary_map.cpp modular_symbol.cpp manin_element.cpp newspace.cpp linalg.cpp)

add_executable(main main.cpp ${SOURCES})
add_library(modular_symbols STATIC ${SOURCES})