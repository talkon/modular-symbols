#ifndef MANIN_BASIS_H
#define MANIN_BASIS_H

// This file defines ManinElement and bases of spaces of Manin symbols

#include "manin_symbol.h"
#include "cache_decorator.h"
#include <flint/fmpq.h>

// Manin generator with coefficient, used as a component in ManinElement.
//
// The generator is stored as the index (in the std::vector returned by ManinGenerator).
// The level N is needed to get an actual ManinGenerator out of an MGWC.
//
// For now, the coefficient is stored by value to simplify memory management.
struct MGWC {
  int64_t index;
  fmpq coeff;

  // Prints this MGWC
  void print();
};

// An element of the space of Manin symbols of level N.
// Represented as a sparse linear combination of basis elements.
struct ManinElement {
  int64_t N;
  std::vector<MGWC> components;

  // Prints this Manin element
  void print();

  // Prints this Manin element, with index expanded into ManinGenerators
  void print_with_generators();
};

// Given a level and an index of a Manin generator
// (within the std:vector returned by `manin_generators`),
// returns that Manin generator as a ManinElement (i.e. a linear combination)
// Results are cached.
ManinElement level_and_index_to_basis(const int64_t level, const int64_t index);

// Computes the basis of the space of Manin symbols of level N.
// Results are cached.
std::vector<ManinGenerator> manin_basis(const int64_t index);

#endif // MANIN_BASIS_H