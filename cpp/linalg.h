#ifndef LINALG_H
#define LINALG_H

// This file defines linear algebra helper functions

#include "manin_basis.h"
#include "manin_element.h"

#include <vector>

// Computes the kernel (represented by a basis) of a given linear map `f` acting on
// a vector space of ManinElements (also represented by a basis).
// `f` should be a map to Manin symbols of level `M`.
std::vector<ManinElement> map_kernel(std::vector<ManinElement>, std::function<ManinElement(ManinBasisElement)> f, int64_t M);

struct DecomposeResult {
  std::vector<std::vector<ManinElement>> done;
  std::vector<std::vector<ManinElement>> remaining;

  static DecomposeResult empty();
};

// Decomposes a subspace with the given basis into simple f-modules,
// where `f` is a linear map.
// If is_atkin_lehner, assume that minimal polynomial divides T^2-1.
DecomposeResult decompose(std::vector<ManinElement>, std::function<ManinElement(ManinBasisElement)> f, bool is_atkin_lehner);


#endif // LINALG_H