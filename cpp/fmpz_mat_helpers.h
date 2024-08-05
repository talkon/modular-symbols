#ifndef FMPZ_MAT_HELPERS_H
#define FMPZ_MAT_HELPERS_H

#include <flint/fmpz_types.h>
#include <flint/fmpq_types.h>

// Computes C = A * B by reducing and multiplying one prime at a time.
void fmpz_mat_mul_nmod(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);

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

// Gets the maximum element of a matrix
void fmpz_mat_max_elt(fmpz_t res, const fmpz_mat_t mat);

// Gets the minimum element of a matrix
void fmpz_mat_min_elt(fmpz_t res, const fmpz_mat_t mat);

// Returns the total amount of memory used by `num` in bytes
ulong fmpz_alloc_size(const fmpz_t num);

#endif