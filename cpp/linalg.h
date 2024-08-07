#ifndef LINALG_H
#define LINALG_H

// This file defines linear algebra helper functions

#include "manin_basis.h"
#include "manin_element.h"
#include "flint_wrappers.h"
#include "subspace_basis.h"

#include <vector>

// Forward declaration
struct Subspace;

// Computes the kernel (represented by a basis) of a given linear map `f` acting on
// a vector space of ManinElements (also represented by a basis).
// `f` should be a map to Manin symbols of level `M`.
DenseBasis map_kernel(DenseBasis&, std::function<ManinElement(ManinBasisElement)> f, int64_t N, int64_t M);

struct SplitResult {
  DenseBasis pos_space;
  DenseBasis neg_space;

  static SplitResult empty(int64_t level);
};

// Splits a space by Atkin-Lehner sign, where `f` is an Atkin-Lehner involution.
SplitResult split(DenseBasis&, std::function<ManinElement (ManinBasisElement)> f, int64_t N);

struct DecomposeResult {
  std::vector<Subspace> done;
  std::vector<Subspace> special;
  std::vector<Subspace> remaining;

  static DecomposeResult empty();
};

// Decomposes a subspace with the given basis into simple f-modules,
// where `map_of_basis` is a matrix of the linear map f acting on the standard basis.
// If `dimension_only` is true, only computes dimension of subspaces.
DecomposeResult decompose(Subspace, FmpqMatrix& map_of_basis, bool dimension_only, bool prime_opt, const slong mem_threshold);

#endif // LINALG_H