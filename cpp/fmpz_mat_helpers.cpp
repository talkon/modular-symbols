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

// Horner's method
void fmpz_poly_apply_fmpq_mat_horner(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f) {
  int degree = fmpz_poly_degree(f);
  assert(degree >= 1);

  fmpz_t coeff;
  fmpz_init(coeff);

  fmpz_poly_get_coeff_fmpz(coeff, f, degree);
  fmpq_mat_scalar_mul_fmpz(dst, src, coeff);

  for (int i = 1; i <= degree; i++) {
    fmpz_poly_get_coeff_fmpz(coeff, f, degree - i);
    for (int j = 0; j < fmpq_mat_nrows(dst); j++) {
      fmpq_add_fmpz(fmpq_mat_entry(dst, j, j), fmpq_mat_entry(dst, j, j), coeff);
    }

    if (i < degree) {
      fmpq_mat_mul(dst, dst, src);
    }
  }

  fmpz_clear(coeff);
}

// Algorithm B in Paterson-Stockmeyer: 2 * sqrt(deg(P)) matrix multiplications
// TODO: consider scaling up everything and work on integer matrices
void fmpz_poly_apply_fmpq_mat_ps(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f) {

  DEBUG_INFO(5,
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

void fmpz_poly_apply_fmpq_mat(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f) {
  int degree = fmpz_poly_degree(f);

  if (degree == 0) {
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_poly_get_coeff_fmpz(coeff, f, 0);
    fmpq_mat_zero(dst);
    for (int j = 0; j < fmpq_mat_nrows(dst); j++) {
      fmpq_set_fmpz(fmpq_mat_entry(dst, j, j), coeff);
    }
    fmpz_clear(coeff);
    return;
  } else if (degree <= 3) {
    fmpz_poly_apply_fmpq_mat_horner(dst, src, f);
    return;
  } else {
    fmpz_poly_apply_fmpq_mat_ps(dst, src, f);
    return;
  }
}

void fmpz_mat_print_dimensions(const fmpz_mat_t mat) {
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

// Code modified from fmpz_mat/max_bits.c in FLINT
slong fmpq_mat_max_bits(const fmpq_mat_t mat) {
  slong i;
  slong bits, row_bits, sign;

  sign = 1;
  bits = 0;

  if (mat->r == 0 || mat->c == 0)
    return 0;

  for (i = 0; i < mat->r; i++)
  {
    row_bits = _fmpz_vec_max_bits((fmpz*) mat->rows[i], 2 * mat->c);
    if (row_bits < 0)
    {
      row_bits = -row_bits;
      sign = -1;
    }
    bits = FLINT_MAX(bits, row_bits);
  }

  return bits * sign;
}

void fmpq_mat_print_dimensions(const fmpq_mat_t mat) {
  int r = fmpq_mat_nrows(mat);
  int c = fmpq_mat_ncols(mat);
  int max = fmpq_mat_max_bits(mat);

  printf("(%d x %d, %d)_Q ", r, c, max);
}



