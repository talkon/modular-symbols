#ifndef MANIN_SYMBOL_H
#define MANIN_SYMBOL_H

// This file defines ManinSymbol, ManinGenerator, and related functions

#include <vector>

// Forward declarations
struct ManinGenerator;
struct ManinElement;
struct ModularSymbol;
struct IntMatrix2x2;

// Represents a Manin symbol. The level N should fit in 32 bits.
struct ManinSymbol {
  int64_t c;
  int64_t d;
  int64_t N;

  // friend auto operator<=> (const ManinSymbol&, const ManinSymbol&) = default;
  bool operator<(const ManinSymbol&) const;

  // Prints this Manin symbol
  void print() const;

  // Checks equivalence between two Manin symbols.
  // Note that this is different from equality (==), which is used for caching.
  bool is_equivalent(const ManinSymbol&) const;

  // Returns (c mod N, d mod N)
  ManinSymbol repr();

  // Applies the eta relation to this Manin symbol. [Cremona Ch 2.5]
  ManinSymbol apply_eta();

  // Applies the S relation to this Manin symbol. [Stein Ch 3.3]
  ManinSymbol apply_S();

  // Applies the T relation to this Manin symbol.  [Stein Ch 3.3]
  ManinSymbol apply_T();

  // Applies the T^2 relation to this Manin symbol.  [Stein Ch 3.3]
  ManinSymbol apply_T_2();

  // Applies the right action of a 2x2 matrix on this Manin symbol. [Stein Ch 3.3]
  ManinSymbol right_action_by(IntMatrix2x2);

  // Converts this Manin symbol to a modular symbol.
  ModularSymbol as_modular_symbol();

  // Returns a ManinGenerator equivalent to this Manin symbol.
  ManinGenerator as_generator();
};

// Represents a Manin generator, i.e. a canonical Manin symbol
// for each element of P^1(Z/NZ)
struct ManinGenerator : ManinSymbol {
  // Index of the Manin generator, in the order provided by `manin_generators()`.
  int64_t index = -1;

  explicit ManinGenerator (int64_t index, ManinSymbol ms): ManinSymbol(ms), index(index) {};
};

// Computes the Manin generators of a given level.
// Returns the generators ordered by (c, d mod N/c).
// Results are cached.
std::vector<ManinGenerator>& manin_generators(const int64_t level);

// Finds the Manin generator that is equivalent to the given Manin symbol.
// Results are cached.
ManinGenerator find_generator(const ManinSymbol);

#endif // MANIN_SYMBOL_H
