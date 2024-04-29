#ifndef HEILBRONN_H
#define HEILBRONN_H

#include "modular_symbol.h"

#include <vector>

// Computes Cremona's Heilbronn matrices for a given prime p.
// See Cremona Ch 2.4
std::vector<IntMatrix2x2>& heilbronn_cremona(int64_t p);

// Computes Merel's Heilbronn matrices for a given prime p.
// See p.87 in Merel
std::vector<IntMatrix2x2>& heilbronn_merel(int64_t p);

#endif // HEILBRONN_H