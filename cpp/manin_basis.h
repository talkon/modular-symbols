#ifndef MANIN_BASIS_H
#define MANIN_BASIS_H

// This file defines ManinElement and bases of spaces of Manin symbols

#include "manin_symbol.h"
#include "cache_decorator.h"
#include <flint/fmpq.h>

// Given a level and an index of a Manin generator
// (within the std:vector returned by `manin_generators`),
// returns that Manin generator as a ManinElement (i.e. a linear combination)
// Results are cached.
ManinElement level_and_index_to_basis(const int64_t level, const int64_t index);

// Computes the basis of the space of Manin symbols of level N.
// Results are cached.
std::vector<ManinGenerator> manin_basis(const int64_t index);

// Converts a modular symbol {a/b} to a ManinElement with the given level.
ManinElement fraction_to_manin_element(const int64_t a, const int64_t b, const int64_t level);

#endif // MANIN_BASIS_H