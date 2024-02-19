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

  friend auto operator<=> (const MGWC& left, const MGWC& right) {
    return left.index <=> right.index;
  }

  // Prints this MGWC
  void print();
};

// An element of the space of Manin symbols of level N.
// Represented as a sparse linear combination of basis elements.
//
// The indices in `components` must be distinct, and if `is_sorted` is true,
// the indices must be sorted.
//
// XXX: this representation seems inefficient for adding ManinElements;
// perhaps consider other representations? (like a map from index -> coeff)
struct ManinElement {
  int64_t N;
  std::vector<MGWC> components;
  bool is_sorted = false;

  // Zero element of a given level
  static ManinElement zero(const int64_t level);

  // Sorts the components of this element.
  void sort();

  // Marks the components of this element as sorted, without checking.
  void mark_as_sorted_unchecked();

  // Overloaded arithmetic operations

  ManinElement& operator+= (const ManinElement&);
  ManinElement& operator-= (const ManinElement&);
  friend ManinElement operator+ (const ManinElement&, const ManinElement&);
  friend ManinElement operator- (const ManinElement&, const ManinElement&);

  // Prints this Manin element
  void print() const;

  // Prints this Manin element, with index expanded into ManinGenerators
  void print_with_generators() const;
};

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