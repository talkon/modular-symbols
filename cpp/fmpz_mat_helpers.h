#ifndef FMPZ_MAT_HELPERS_H
#define FMPZ_MAT_HELPERS_H

#include <flint/fmpz_types.h>
#include <flint/fmpq_types.h>

// Sets dst to f(src).
// `dst` and `src` must be square matrices with the same dimension, and cannot alias.
void fmpq_poly_apply_fmpq_mat(fmpq_mat_t dst, const fmpq_mat_t src, const fmpq_poly_t f);

// Sets dst to f(src) using Horner's method.
// `dst` and `src` must be square matrices with the same dimension, and cannot alias.
void fmpz_poly_apply_fmpq_mat_horner(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f);

// Sets dst to f(src) using Paterson-Stockmeyer's (non-recursive) algorithm.
// `dst` and `src` must be square matrices with the same dimension, and cannot alias.
void fmpz_poly_apply_fmpq_mat_ps(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f);

// Divides each row by the gcd of the row.
void fmpz_mat_div_rowwise_gcd(fmpz_mat_t mat);

// Divides each row by the gcd of the row.
void fmpz_mat_div_colwise_gcd(fmpz_mat_t mat);

// Computes nullspace using fmpz_mat_rref_mul
slong fmpz_mat_nullspace_mul(fmpz_mat_t res, const fmpz_mat_t mat);

// Prints (row x col, max_bits) of a matrix.
void fmpz_mat_print_dimensions(const fmpz_mat_t mat);

// Prints (row x col, max_bits) of a matrix.
void fmpq_mat_print_dimensions(const fmpq_mat_t mat);

#endif