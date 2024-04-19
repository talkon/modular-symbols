#ifndef HEILBRONN_H
#define HEILBRONN_H

#include "modular_symbol.h"

#include <vector>

// Computes Cremona's Heilbronn matrices for a given prime p.
std::vector<IntMatrix2x2> heilbronn_cremona(int64_t p);

#endif // HEILBRONN_H