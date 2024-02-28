#ifndef NEWFORM_SUBSPACES_H
#define NEWFORM_SUBSPACES_H

#include "modular_symbol.h"

#include <vector>

// Forward declarations
struct ManinBasisElement;
struct ManinElement;

// Computes the Heilbronn matrices for a given prime p.
std::vector<IntMatrix2x2> heilbronn_matrices(int64_t p);

// Computes the action of the Hecke operator T_p on the given Manin basis element.
ManinElement hecke_action(ManinBasisElement, int64_t p);

// Computes the newform subspaces of a given level.
std::vector<std::vector<ManinElement>> newform_subspaces(int64_t level);

// Computes the dimensions of the newform subspaces of a given level.
std::vector<int> newform_subspace_dimensions(int64_t level);

#endif // NEWFORM_SUBSPACES_H