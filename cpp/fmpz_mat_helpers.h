#ifndef FMPZ_MAT_HELPERS_H
#define FMPZ_MAT_HELPERS_H

#include <flint/fmpz_types.h>
#include <flint/fmpq_types.h>

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

// Returns the total amount of memory used by `num` in bytes
ulong fmpz_alloc_size(const fmpz_t num);

#endif