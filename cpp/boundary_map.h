#ifndef BOUNDARY_MAP_H
#define BOUNDARY_MAP_H

// This file defines the boundary map, as in [Stein Ch 3.5]

#include "manin_basis.h"
#include "manin_element.h"
#include "utils.h"

// Represents a cusp {c/d} of level N
struct Cusp {
  int64_t c;
  int64_t d;
  int64_t N;

  friend auto operator<=> (const Cusp&, const Cusp&) = default;

  // Prints this cusp.
  void print();

  // Computes xgcd(c, d)
  utils::XgcdResult xgcd() const;

  // Returns the negation {-c/d} of the cusp
  Cusp negated() const;

  // Checks equivalence between two cusps.
  // Note that this is different from equality (==), which is used for caching.
  bool is_equivalent(const Cusp&);
};

// Computes the basis of the kernel of the boundary map for Manin symbols of level N,
// i.e. the basis of the space of cuspidal modular symbols of level N.
std::vector<ManinElement> cuspidal_manin_basis(int64_t);

#endif // BOUNDARY_MAP_H