#ifndef MATRIX_POLY_H
#define MATRIX_POLY_H

#include <flint/fmpz_types.h>
#include <flint/fmpq_types.h>

// Sets dst to f(src) using Horner's method.
// `dst` and `src` must be square matrices with the same dimension, and cannot alias.
void fmpz_poly_apply_fmpq_mat_horner(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f);

// Sets dst to f(src) using Paterson-Stockmeyer's (non-recursive) algorithm.
// `dst` and `src` must be square matrices with the same dimension, and cannot alias.
void fmpz_poly_apply_fmpq_mat_ps(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f, const slong mem_threshold);

// Sets dst to f(src), selecting between Horner's and P-S method.
// `dst` and `src` must be square matrices with the same dimension, and cannot alias.
void fmpz_poly_apply_fmpq_mat(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f, const slong mem_threshold);

#endif // MATRIX_POLY_H