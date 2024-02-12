#ifndef MANIN_SYMBOL_H
#define MANIN_SYMBOL_H

#include "cache_decorator.h"
#include <flint/fmpq.h>

// Forward declarations
struct ManinGenerator;
struct ManinElement;

// Represents a Manin symbol. The level N should fit in 32 bits.
struct ManinSymbol {
  int64_t c;
  int64_t d;
  int64_t N;

  friend auto operator<=> (const ManinSymbol&, const ManinSymbol&) = default;

  // Prints this Manin symbol
  void print();

  // Checks equivalence between two Manin symbols.
  // Note that this is different from equality (==), which is used for caching.
  bool is_equivalent(const ManinSymbol&);

  // Applies the eta relation to this Manin symbol. [Cremona Ch 2.5]
  ManinSymbol apply_eta();

  // Applies the S relation to this Manin symbol. [Stein Ch 3.3]
  ManinSymbol apply_S();

  // Applies the T relation to this Manin symbol.  [Stein Ch 3.3]
  ManinSymbol apply_T();

  // Applies the T^2 relation to this Manin symbol.  [Stein Ch 3.3]
  ManinSymbol apply_T_2();

  // Returns a pointer to the generator equivalent to this Manin symbol.
  ManinGenerator as_generator();
};

// Represents a Manin generator, i.e. a canonical Manin symbol
// for each element of P^1(Z/NZ)
struct ManinGenerator : ManinSymbol {
  // Index of the Manin generator, in the order provided by `manin_generators()`.
  int64_t index = -1;

  explicit ManinGenerator (int64_t index, ManinSymbol ms): ManinSymbol(ms), index(index) {};

  // Converts to a Manin element, assuming this generator is in the basis.
  ManinElement as_element_unchecked();
};

// Computes the Manin generators of a given level.
// Returns the generators ordered by (c, d mod N/c).
// Results are cached.
inline std::vector<ManinGenerator> manin_generators(const int64_t level);

// Finds the Manin generator that is equivalent to the given Manin symbol.
// Results are cached.
inline ManinGenerator find_generator(const ManinSymbol);

// Given a level and an index of a Manin generator
// (within the std:vector returned by `manin_generators`),
// returns that Manin generator as a ManinElement (i.e. a linear combination)
inline ManinElement level_and_index_to_basis(const int64_t level, const int64_t index);

// Computes the basis of the space of Manin symbols of level N.
inline std::vector<ManinGenerator> manin_basis(const int64_t index);

// Manin generator with coefficient, used as a component in ManinElement.

// The generator is stored as the index (in the std::vector returned by ManinGenerator).
// The level N is needed to get an actual ManinGenerator out of an MGWC.

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

// typedef struct {
//   int64_t a;
//   int64_t b;
//   int64_t c;
//   int64_t d;
// } modular_symbol_t;

// typedef struct {
//   int64_t a;
//   int64_t b;
//   int64_t c;
//   int64_t d;
// } sl2z_element_t;

#endif // MANIN_SYMBOL_H
