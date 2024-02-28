#ifndef LINALG_H
#define LINALG_H

// This file defines linear algebra helper functions

#include "manin_basis.h"
#include "manin_element.h"

#include <flint/fmpz_types.h>
#include <flint/fmpq_types.h>

#include <vector>

// Computes the kernel (represented by a basis) of a given linear map `f` acting on
// a vector space of ManinElements (also represented by a basis).
// `f` should be a map to Manin symbols of level `M`.
std::vector<ManinElement> map_kernel(std::vector<ManinElement>, std::function<ManinElement(ManinBasisElement)> f, int64_t M);

// Sets dst to f(src).
// `dst` and `src` must be square matrices with the same dimension, and cannot alias.
void fmpq_poly_apply_fmpq_mat(fmpq_mat_t dst, const fmpq_mat_t src, const fmpq_poly_t f);

// Sets dst to f(src).
// `dst` and `src` must be square matrices with the same dimension, and cannot alias.
void fmpz_poly_apply_fmpq_mat(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f);

struct DecomposeResult {
  std::vector<std::vector<ManinElement>> done;
  std::vector<std::vector<ManinElement>> remaining;

  static DecomposeResult empty();
};

// Decomposes a subspace with the given basis into simple f-modules,
// where `f` is a linear map.
DecomposeResult decompose(std::vector<ManinElement>, std::function<ManinElement(ManinBasisElement)> f);


#endif // LINALG_H