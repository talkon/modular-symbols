#include "fmpz_mat_helpers.h"
#include "debug_utils.h"

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/ulong_extras.h>

#include <array>
#include <cassert>

void fmpz_mat_div_rowwise_gcd(fmpz_mat_t mat) {
  fmpz_mat_t window;
  fmpz_t den;
  fmpz_init(den);

  int nrows = fmpz_mat_nrows(mat);
  int ncols = fmpz_mat_ncols(mat);

  for (int row = 0; row < nrows; row++) {
    fmpz_mat_window_init(window, mat, row, 0, row+1, ncols);

    fmpz_mat_content(den, window);
    if (!fmpz_is_zero(den)) {
      fmpz_mat_scalar_divexact_fmpz(window, window, den);
    }
    fmpz_mat_window_clear(window);
  }

  fmpz_clear(den);
}

void fmpz_mat_div_colwise_gcd(fmpz_mat_t mat) {
  fmpz_mat_t window;
  fmpz_t den;
  fmpz_init(den);

  int nrows = fmpz_mat_nrows(mat);
  int ncols = fmpz_mat_ncols(mat);

  for (int col = 0; col < ncols; col++) {
    fmpz_mat_window_init(window, mat, 0, col, nrows, col+1);

    fmpz_mat_content(den, window);
    if (!fmpz_is_zero(den)) {
      fmpz_mat_scalar_divexact_fmpz(window, window, den);
    }
    fmpz_mat_window_clear(window);
  }

  fmpz_clear(den);
}

void fmpq_poly_apply_fmpq_mat(fmpq_mat_t dst, const fmpq_mat_t src, const fmpq_poly_t f) {
  fmpq_t coeff;
  fmpq_init(coeff);
  int degree = fmpq_poly_degree(f);

  fmpq_mat_zero(dst);

  // XXX: might be faster to just add to diagonal elements directly.
  fmpq_mat_t one, coeff_times_one;
  fmpq_mat_init_set(one, dst);
  fmpq_mat_one(one);
  fmpq_mat_init_set(coeff_times_one, one);

  for (int i = 0; i <= degree; i++) {
    fmpq_poly_get_coeff_fmpq(coeff, f, degree - i);
    fmpq_mat_scalar_mul_fmpq(coeff_times_one, one, coeff);
    fmpq_mat_add(dst, dst, coeff_times_one);
    if (i != degree) {
      fmpq_mat_mul(dst, dst, src);
    }
  }

  fmpq_mat_clear(one);
  fmpq_mat_clear(coeff_times_one);
  fmpq_clear(coeff);
}

// Horner's method
void fmpz_poly_apply_fmpq_mat_horner(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f) {
  fmpz_t coeff;
  fmpz_init(coeff);
  int degree = fmpz_poly_degree(f);

  DEBUG_INFO(3,
    {
      printf("fmpz_poly_apply_fmpq_mat_horner called with f(T) = ");
      fmpz_poly_print_pretty(f, "T");
      printf("\n");
    }
  )

  fmpq_mat_zero(dst);

  // XXX: might be faster to just add to diagonal elements directly.
  fmpq_mat_t one, coeff_times_one;
  fmpq_mat_init_set(one, dst);
  fmpq_mat_one(one);
  fmpq_mat_init_set(coeff_times_one, one);

  for (int i = 0; i <= degree; i++) {
    fmpz_poly_get_coeff_fmpz(coeff, f, degree - i);
    fmpq_mat_scalar_mul_fmpz(coeff_times_one, one, coeff);
    fmpq_mat_add(dst, dst, coeff_times_one);
    if (i != degree) {
      fmpq_mat_mul(dst, dst, src);
    }
  }

  fmpq_mat_clear(one);
  fmpq_mat_clear(coeff_times_one);
  fmpz_clear(coeff);
}

// Algorithm B in Paterson-Stockmeyer: 2 * sqrt(deg(P)) matrix multiplications
// TODO: consider scaling up everything and work on integer matrices
void fmpz_poly_apply_fmpq_mat_ps(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f) {

  DEBUG_INFO(3,
    {
      printf("fmpz_poly_apply_fmpq_mat_ps called with f(T) = ");
      fmpz_poly_print_pretty(f, "T");
      printf("\n");
    }
  )

  ulong l = fmpz_poly_degree(f);
  ulong k = n_sqrt(l);

  ulong d = fmpq_mat_nrows(src);
  assert(d == fmpq_mat_ncols(src));

  // Compute T, T^2, .., T^k
  fmpq_mat_t* pows = (fmpq_mat_t*) flint_malloc((k + 1) * sizeof(fmpq_mat_t));
  fmpq_mat_t tmp;

  fmpq_mat_init(pows[0], d, d);
  fmpq_mat_one(pows[0]);

  fmpq_mat_init_set(tmp, src);
  for (int i = 1; i < k; i++) {
    fmpq_mat_init_set(pows[i], tmp);
    fmpq_mat_mul(tmp, tmp, src);
  }

  fmpq_mat_init_set(pows[k], tmp);

  fmpq_mat_init(dst, d, d);
  fmpq_mat_zero(dst);

  fmpz_t a;

  fmpz_init(a);
  for (int j = l / k; j >= 0; j--) {

    for (int i = 0; i < k; i++) {
      fmpz_poly_get_coeff_fmpz(a, f, j * k + i);
      if (!fmpz_is_zero(a)) {
        fmpq_mat_scalar_mul_fmpz(tmp, pows[i], a);
        fmpq_mat_add(dst, dst, tmp);
      }
    }

    if (j > 0) {
      fmpq_mat_mul(dst, dst, pows[k]);
    }
  }

  fmpz_clear(a);
  fmpq_mat_clear(tmp);
  for (int i = 0; i <= k; i++) {
    fmpq_mat_clear(pows[i]);
  }
  flint_free(pows);
}

void fmpz_mat_print_dimensions(fmpz_mat_t mat) {
  int r = fmpz_mat_nrows(mat);
  int c = fmpz_mat_ncols(mat);
  int max = fmpz_mat_max_bits(mat);

  printf("(%d x %d, %d)", r, c, max);
}

// Code taken from fmpz_mat/nullspace.c in FLINT, modified to use fmpz_mat_rref_mul
slong fmpz_mat_nullspace_mul(fmpz_mat_t res, const fmpz_mat_t mat) {
  slong i, j, k, n, rank, nullity;
  slong * pivots;
  slong * nonpivots;
  fmpz_mat_t tmp;
  fmpz_t den;

  n = mat->c;

  fmpz_mat_init_set(tmp, mat);
  fmpz_init(den);

  rank = fmpz_mat_rref_mul(tmp, den, mat);
  nullity = n - rank;

  fmpz_mat_zero(res);
  if (rank == 0) {
    for (i = 0; i < nullity; i++)
      fmpz_one(res->rows[i] + i);
  }
  else if (nullity) {
    pivots = (slong*) flint_malloc(rank * sizeof(slong));
    nonpivots = (slong*) flint_malloc(nullity * sizeof(slong));

    for (i = j = k = 0; i < rank; i++) {
      while (fmpz_is_zero(tmp->rows[i] + j)) {
        nonpivots[k] = j;
        k++;
        j++;
      }
      pivots[i] = j;
      j++;
    }
    while (k < nullity) {
      nonpivots[k] = j;
      k++;
      j++;
    }

    fmpz_set(den, tmp->rows[0] + pivots[0]);

    for (i = 0; i < nullity; i++) {
      for (j = 0; j < rank; j++)
        fmpz_set(res->rows[pivots[j]] + i, tmp->rows[j] + nonpivots[i]);
      fmpz_neg(res->rows[nonpivots[i]] + i, den);
    }

    flint_free(pivots);
    flint_free(nonpivots);
  }

  fmpz_clear(den);
  fmpz_mat_clear(tmp);

  return nullity;
}
