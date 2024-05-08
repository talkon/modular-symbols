set(SOURCES manin_symbol.cpp manin_basis.cpp boundary_map.cpp modular_symbol.cpp manin_element.cpp newspace.cpp linalg.cpp newform_subspaces.cpp debug_utils.cpp fmpz_mat_helpers.cpp heilbronn.cpp flint_wrappers.cpp subspace.cpp)

add_executable(main main.cpp ${SOURCES})
add_library(modular_symbols STATIC ${SOURCES})