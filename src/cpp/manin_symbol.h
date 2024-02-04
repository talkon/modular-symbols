#ifndef MANIN_SYMBOL_H
#define MANIN_SYMBOL_H

#include "cache_decorator.h"

/**
 * Represents a Manin symbol. The level N should fit in 32 bits.
 */
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
};

struct ManinGenerator : ManinSymbol {
  explicit ManinGenerator (ManinSymbol ms): ManinSymbol(ms) {};
};

// Computes the Manin generators of level N.
// Returns the generators ordered by (c, d mod N/c).
// Results are cached.
inline std::vector<ManinGenerator> manin_generators(const int64_t);

// Finds the index of the Manin generator (in the order provided by `manin_generators()`)
// that is equivalent to the given Manin symbol.
// Results are cached.
inline int64_t find_generator_index(const ManinSymbol&);

// Finds the Manin generator that is equivalent to the given Manin symbol.
inline ManinGenerator find_generator(const ManinSymbol&);

// /**
//  * Manin generator with coefficient
// */
// struct MGWC {
//   ManinGenerator generator;
//   fmpq_t coeff;
// }; // manin generator with coefficient

// typedef struct {
//   int64_t N;
//   std::vector<MGWC> components;
// } manin_element_t;

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
