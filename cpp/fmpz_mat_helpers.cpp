#include "fmpz_mat_helpers.h"
#include "debug_utils.h"

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/ulong_extras.h>

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

ulong fmpz_alloc_size(const fmpz_t num) {
  if (COEFF_IS_MPZ(*num)) {
    __mpz_struct *z = COEFF_TO_PTR(*num);
    return sizeof(__mpz_struct) + z->_mp_alloc * 8;
  }
  else return 8;
}

ulong fmpz_mat_total_size(const fmpz_mat_t mat) {
  ulong size = 0;
  for (slong i = 0; i < mat->r; i++) {
    for (slong j = 0; j < mat->c; j++) {
      size += fmpz_alloc_size(((fmpz*) mat->rows[i]) + j);
    }
  }
  return size + sizeof(fmpz_mat_struct);
}

ulong fmpq_mat_total_size(const fmpq_mat_t mat) {
  ulong size = 0;
  for (slong i = 0; i < mat->r; i++) {
    for (slong j = 0; j < 2 * mat->c; j++) {
      size += fmpz_alloc_size(((fmpz*) mat->rows[i]) + j);
    }
  }
  return size + sizeof(fmpq_mat_struct);
}


void fmpz_mat_print_dimensions(const fmpz_mat_t mat) {
  int r = fmpz_mat_nrows(mat);
  int c = fmpz_mat_ncols(mat);
  int max = fmpz_mat_max_bits(mat);
  ulong size = fmpz_mat_total_size(mat);

  printf("(%d x %d, %d, %lu)", r, c, max, size);
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

void fmpz_mat_max_elt(fmpz_t res, const fmpz_mat_t mat) {
  assert(mat->r > 0 && mat->c > 0);

  fmpz_set(res, fmpz_mat_entry(mat, 0, 0));

  for (slong i = 0; i < mat->r; i++) {
    for (slong j = 0; j < mat->c; j++) {
      if (fmpz_cmp(res, fmpz_mat_entry(mat, i, j)) > 0) {
        fmpz_set(res, fmpz_mat_entry(mat, i, j));
      }
    }
  }
}

void fmpz_mat_min_elt(fmpz_t res, const fmpz_mat_t mat) {
  assert(mat->r > 0 && mat->c > 0);

  fmpz_set(res, fmpz_mat_entry(mat, 0, 0));

  for (slong i = 0; i < mat->r; i++) {
    for (slong j = 0; j < mat->c; j++) {
      if (fmpz_cmp(res, fmpz_mat_entry(mat, i, j)) < 0) {
        fmpz_set(res, fmpz_mat_entry(mat, i, j));
      }
    }
  }
}

void fmpq_mat_print_dimensions(const fmpq_mat_t mat) {
  int r = fmpq_mat_nrows(mat);
  int c = fmpq_mat_ncols(mat);
  int max = fmpq_mat_max_bits(mat);
  ulong size = fmpq_mat_total_size(mat);

  printf("(%d x %d, %d, %lu)_Q ", r, c, max, size);
}



