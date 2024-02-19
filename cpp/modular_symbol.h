#ifndef MODULAR_SYMBOL_H
#define MODULAR_SYMBOL_H

// This file defines ModularSymbol and related functions.

#include <vector>

// Forward declaration
struct ManinElement;

// Represents a matrix ((x y) (z w)) of integers.
struct IntMatrix2x2 {
  int64_t x;
  int64_t y;
  int64_t z;
  int64_t w;
};

// Represents the modular symbol {a/b, c/d}.
struct ModularSymbol {
  int64_t a;
  int64_t b;
  int64_t c;
  int64_t d;

  ModularSymbol left_action_by(const IntMatrix2x2);

  // Converts this modular symbol to a ManinElement of a given level.
  ManinElement to_manin_element(const int64_t level);
};



#endif // MODULAR_SYMBOL_H