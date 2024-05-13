#ifndef NEWFORM_SUBSPACES_H
#define NEWFORM_SUBSPACES_H

#include "modular_symbol.h"
#include "manin_element.h"
#include "flint_wrappers.h"
#include "subspace.h"

#include <vector>

// Forward declarations
struct ManinBasisElement;

// Computes the action of the Hecke operator T_p on the given Manin basis element.
// p must not divide the level N.
ManinElement& hecke_action(ManinBasisElement, int64_t p);

// Computes the matrix of the action of the Hecke operator T_p.
// p must not divide the level.
FmpqMatrix hecke_matrix(int64_t level, int64_t p);

// Computes the action of the Atkin-Lehner involution w_q on the given Manin basis element.
// q must be a prime power such that q || N.
ManinElement& atkin_lehner_action(ManinBasisElement, int64_t q);

// Computes the newform subspaces of a given level.
std::vector<Subspace> newform_subspaces(int64_t level, bool dimension_only, int trace_depth, bool prime_opt);

// Computes the dimensions of the newform subspaces of a given level.
std::vector<int> newform_subspace_dimensions(int64_t level);

#endif // NEWFORM_SUBSPACES_H