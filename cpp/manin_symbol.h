#ifndef MANIN_SYMBOL_H
#define MANIN_SYMBOL_H

// This file defines ManinSymbol, ManinGenerator, and related functions

#include <vector>

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
std::vector<ManinGenerator> manin_generators(const int64_t level);

// Finds the Manin generator that is equivalent to the given Manin symbol.
// Results are cached.
ManinGenerator find_generator(const ManinSymbol);

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
