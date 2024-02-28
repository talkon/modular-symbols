#ifndef NEWFORM_SUBSPACES_H
#define NEWFORM_SUBSPACES_H

#include "modular_symbol.h"

// Computes the Heilbronn matrices for a given prime p.
std::vector<IntMatrix2x2> heilbronn_matrices(int64_t p);

#endif // NEWFORM_SUBSPACES_H