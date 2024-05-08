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

  // Compute 1, T, T^2, .., T^k
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

void fmpz_poly_apply_fmpz_mat_ps(fmpz_mat_t dst, const fmpz_mat_t src, const fmpz_poly_t f) {

  DEBUG_INFO(5,
    {
      printf("fmpz_poly_apply_fmpq_mat_ps called with f(T) = ");
      fmpz_poly_print_pretty(f, "T");
      printf("\n");
    }
  )

  ulong l = fmpz_poly_degree(f);
  ulong k = n_sqrt(l);

  ulong d = fmpz_mat_nrows(src);
  assert(d == fmpz_mat_ncols(src));

  // Compute 1, T, T^2, .., T^k
  fmpz_mat_t* pows = (fmpz_mat_t*) flint_malloc((k + 1) * sizeof(fmpz_mat_t));
  fmpz_mat_t tmp;

  fmpz_mat_init(pows[0], d, d);
  fmpz_mat_one(pows[0]);

  fmpz_mat_init_set(tmp, src);
  for (int i = 1; i < k; i++) {
    fmpz_mat_init_set(pows[i], tmp);
    fmpz_mat_mul(tmp, tmp, src);
  }

  fmpz_mat_init_set(pows[k], tmp);

  fmpz_mat_init(dst, d, d);
  fmpz_mat_zero(dst);

  fmpz_t a;

  fmpz_init(a);
  for (int j = l / k; j >= 0; j--) {

    for (int i = 0; i < k; i++) {
      fmpz_poly_get_coeff_fmpz(a, f, j * k + i);
      if (!fmpz_is_zero(a)) {
        fmpz_mat_scalar_mul_fmpz(tmp, pows[i], a);
        fmpz_mat_add(dst, dst, tmp);
      }
    }

    if (j > 0) {
      fmpz_mat_mul(dst, dst, pows[k]);
    }
  }

  fmpz_clear(a);
  fmpz_mat_clear(tmp);
  for (int i = 0; i <= k; i++) {
    fmpz_mat_clear(pows[i]);
  }
  flint_free(pows);
}

// Algorithm B, but clearing denominators first and applying the algorithm to a fmpz_mat_t.
// This is slower in practice.
void fmpz_poly_apply_fmpq_mat_ps_clear_denom(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f) {

  ulong d = fmpq_mat_nrows(src);
  assert(d == fmpq_mat_ncols(src));

  fmpz_mat_t src_z, dst_z;
  fmpz_t den;
  fmpz_mat_init(src_z, d, d);
  fmpz_mat_init(dst_z, d, d);
  fmpz_init(den);

  fmpq_mat_get_fmpz_mat_matwise(src_z, den, src);

  fmpz_poly_t f_adj;
  fmpz_poly_init(f_adj);
  fmpz_poly_set(f_adj, f);

  int deg = f->length - 1;
  fmpz_t den_pow;
  fmpz_init_set(den_pow, den);

  for (int i = deg - 1; i >= 0; i--) {
    fmpz_mul(fmpz_poly_get_coeff_ptr(f_adj, i), fmpz_poly_get_coeff_ptr(f_adj, i), den_pow);
    if (i > 0) fmpz_mul(den_pow, den_pow, den);
  }

  fmpz_poly_apply_fmpz_mat_ps(dst_z, src_z, f_adj);
  fmpq_mat_set_fmpz_mat_div_fmpz(dst, dst_z, den_pow);

  fmpz_poly_clear(f_adj);
  fmpz_clear(den);
  fmpz_clear(den_pow);
  fmpz_mat_clear(src_z);
  fmpz_mat_clear(dst_z);
}

void _fmpq_poly_apply_fmpq_mat_base(
  fmpq_mat_t dst,
  const fmpq_mat_t src,
  const fmpq_poly_t f,
  const fmpq_mat_t* pows,
  ulong k
) {

  DEBUG_INFO(5,
    {
      printf("_fmpq_poly_apply_fmpq_mat_base called with f(T) = ");
      fmpq_poly_print_pretty(f, "T");
      printf("\n");
    }
  )

  assert(fmpq_poly_is_canonical(f));
  ulong l = fmpq_poly_degree(f);
  assert(l <= k);

  ulong d = fmpq_mat_nrows(src);
  assert(d == fmpq_mat_ncols(src));

  fmpq_mat_t temp;
  fmpq_mat_init(temp, d, d);

  fmpq_t coeff;
  fmpq_init(coeff);

  fmpq_mat_zero(dst);

  for (int i = 0; i <= l; i++) {
    fmpq_poly_get_coeff_fmpq(coeff, f, i);
    fmpq_mat_scalar_mul_fmpq(temp, *(pows + i), coeff);
    fmpq_mat_add(dst, dst, temp);
  }

  fmpq_mat_clear(temp);
  fmpq_clear(coeff);
}

void _fmpq_poly_apply_fmpq_mat_ps_recur_helper(
  fmpq_mat_t dst,
  const fmpq_mat_t src,
  const fmpq_poly_t f,
  const fmpq_mat_t* pows,
  ulong k,
  const fmpq_mat_t* pows_exp,
  ulong m,
  ulong t
) {

  DEBUG_INFO(5,
    {
      printf("_fmpq_poly_apply_fmpq_mat_ps_recur_helper called with f(T) = ");
      fmpq_poly_print_pretty(f, "T");
      printf("\n");
    }
  )

  assert(fmpq_poly_is_canonical(f));
  assert(t >= 1);

  ulong l = fmpq_poly_degree(f);
  ulong p = 1 << (t - 1);

  ulong d = fmpq_mat_nrows(src);
  assert(d == fmpq_mat_ncols(src));

  assert(t <= m);
  assert(l == k * (2 * p - 1));

  fmpq_t leading_coeff;
  fmpq_init(leading_coeff);
  fmpq_poly_get_coeff_fmpq(leading_coeff, f, l);
  assert(fmpq_is_one(leading_coeff));
  fmpq_clear(leading_coeff);

  if (t == 1) {
    _fmpq_poly_apply_fmpq_mat_base(dst, src, f, pows, k);
    return;
  }

  fmpq_poly_t q, r;
  fmpq_poly_init(q);
  fmpq_poly_init(r);

  fmpq_poly_set_trunc(r, f, k * p);
  fmpq_poly_get_slice(q, f, k * p, k * (2 * p - 1));
  fmpq_poly_shift_right(q, q, k * p);

  fmpq_poly_canonicalise(r);
  fmpq_poly_canonicalise(q);

  fmpq_t temp;
  fmpq_init(temp);
  fmpq_poly_get_coeff_fmpq(temp, r, k * (p - 1));
  fmpq_sub_si(temp, temp, 1);
  fmpq_poly_set_coeff_fmpq(r, k * (p - 1), temp);

  fmpq_poly_t c, s;
  fmpq_poly_init(c);
  fmpq_poly_init(s);

  fmpq_poly_divrem(c, s, r, q);

  fmpq_poly_set_coeff_si(s, k * (p - 1), 1);

  fmpq_mat_t q_mat;
  fmpq_mat_init(q_mat, d, d);

  _fmpq_poly_apply_fmpq_mat_ps_recur_helper(q_mat, src, q, pows, k, pows_exp, m, t - 1);

  _fmpq_poly_apply_fmpq_mat_base(dst, src, c, pows, k);
  fmpq_mat_add(dst, dst, *(pows_exp + t));
  fmpq_mat_mul(dst, dst, q_mat);

  fmpq_mat_clear(q_mat);

  fmpq_mat_t s_mat;
  fmpq_mat_init(s_mat, d, d);

  _fmpq_poly_apply_fmpq_mat_ps_recur_helper(s_mat, src, s, pows, k, pows_exp, m, t - 1);
  fmpq_mat_add(dst, dst, s_mat);

  fmpq_mat_clear(s_mat);

  fmpq_poly_clear(q);
  fmpq_poly_clear(r);

  fmpq_clear(temp);

  fmpq_poly_clear(c);
  fmpq_poly_clear(s);
}

void _fmpq_poly_apply_fmpq_mat_ps_recur(
  fmpq_mat_t dst,
  const fmpq_mat_t src,
  const fmpq_poly_t f,
  const fmpq_mat_t* pows,
  ulong k,
  const fmpq_mat_t* pows_exp,
  ulong m
) {

  DEBUG_INFO(5,
    {
      printf("_fmpq_poly_apply_fmpq_mat_ps_recur called with f(T) = ");
      fmpq_poly_print_pretty(f, "T");
      printf("\n");
    }
  )

  if (fmpq_poly_is_zero(f)) {
    fmpq_mat_zero(dst);
    return;
  }

  ulong d = fmpq_mat_nrows(src);
  assert(d == fmpq_mat_ncols(src));

  ulong l = fmpq_poly_degree(f);
  if (l <= k) {
    _fmpq_poly_apply_fmpq_mat_base(dst, src, f, pows, k);
    return;
  } else if (l <= 2 * k) {
    fmpq_poly_t hi, base;
    fmpq_poly_init(hi);
    fmpq_poly_init(base);

    fmpq_poly_set_trunc(base, f, k);
    fmpq_poly_get_slice(hi, f, k, l + 1);
    fmpq_poly_shift_right(hi, hi, k);

    fmpq_mat_t base_mat;
    fmpq_mat_init(base_mat, d, d);
    _fmpq_poly_apply_fmpq_mat_base(base_mat, src, base, pows, k);

    _fmpq_poly_apply_fmpq_mat_base(dst, src, hi, pows, k);
    fmpq_mat_mul(dst, dst, *(pows + k));
    fmpq_mat_add(dst, dst, base_mat);

    fmpq_poly_clear(hi);
    fmpq_poly_clear(base);
    fmpq_mat_clear(base_mat);
    return;
  }

  ulong t = n_flog(l / k, 2);
  assert(t <= m);
  assert(t >= 1);

  DEBUG_INFO(5,
    {
      printf("k: %ld, m: %ld, l: %ld, t: %ld\n", k, m, l, t);
    }
  )


  // if (l == k * (1 << t)) {

  // }

  // fmpq_t leading_coeff;
  // fmpq_init(leading_coeff);
  // fmpq_poly_get_coeff_fmpq(leading_coeff, f, l);

  // fmpq_poly_t f_monic;
  // fmpq_poly_init(f_monic);
  // fmpq_poly_scalar_div_fmpq(f_monic, f, leading_coeff);

  // if (!fmpq_is_one(leading_coeff)) {
  //   fmpq_mat_scalar_mul_fmpq(dst, dst, leading_coeff);
  // }

  // fmpq_poly_clear(f_monic);
  // fmpq_clear(leading_coeff);

  fmpq_poly_t hi, lo, base;
  fmpq_poly_init(hi);
  fmpq_poly_init(lo);
  fmpq_poly_init(base);

  fmpq_poly_set_trunc(base, f, k);

  fmpq_poly_get_slice(lo, f, k, k * (1 << t));
  fmpq_poly_set_coeff_si(lo, k * (1 << t), 1);
  fmpq_poly_shift_right(lo, lo, k);

  fmpq_poly_get_slice(hi, f, k * (1 << t), l + 1);
  fmpq_poly_shift_right(hi, hi, k * (1 << t));
  fmpq_poly_sub_si(hi, hi, 1);

  fmpq_poly_canonicalise(base);
  fmpq_poly_canonicalise(lo);
  fmpq_poly_canonicalise(hi);

  fmpq_mat_t base_mat;
  fmpq_mat_t lo_mat;
  fmpq_mat_init(base_mat, d, d);
  fmpq_mat_init(lo_mat, d, d);

  _fmpq_poly_apply_fmpq_mat_base(base_mat, src, base, pows, k);
  _fmpq_poly_apply_fmpq_mat_ps_recur_helper(lo_mat, src, lo, pows, k, pows_exp, m, t);

  fmpq_mat_mul(lo_mat, lo_mat, *(pows_exp));
  fmpq_mat_add(lo_mat, lo_mat, base_mat);
  fmpq_mat_clear(base_mat);

  _fmpq_poly_apply_fmpq_mat_ps_recur(dst, src, hi, pows, k, pows_exp, m);

  fmpq_mat_mul(dst, dst, *(pows_exp + t));
  fmpq_mat_add(dst, dst, lo_mat);
  fmpq_mat_clear(lo_mat);

  fmpq_poly_clear(hi);
  fmpq_poly_clear(lo);
  fmpq_poly_clear(base);
}

// BUG: using this function seems to give wrong results in some cases, and is also just very slow!
void fmpz_poly_apply_fmpq_mat_ps_recur(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f) {

  DEBUG_INFO(5,
    {
      printf("fmpz_poly_apply_fmpq_mat_ps_recur called with f(T) = ");
      fmpz_poly_print_pretty(f, "T");
      printf("\n");
    }
  )

  ulong l = fmpz_poly_degree(f);

  ulong k = n_sqrt(2 * l);
  ulong m = n_flog(l / k, 2);

  ulong d = fmpq_mat_nrows(src);
  assert(d == fmpq_mat_ncols(src));

  // Compute 1, T, T^2, .., T^k, T^{2k}, ..., T^{k2^m}

  // For 0 <= i <= k, pows[i] = T^i
  fmpq_mat_t* pows = (fmpq_mat_t*) flint_malloc((k + m + 1) * sizeof(fmpq_mat_t));

  // For 0 <= i <= m, pows_exp[i] = T^{2^m k}
  fmpq_mat_t* pows_exp = pows + k;

  fmpq_mat_t tmp;

  fmpq_mat_init(pows[0], d, d);
  fmpq_mat_one(pows[0]);

  fmpq_mat_init_set(tmp, src);
  for (int i = 1; i < k; i++) {
    fmpq_mat_init_set(pows[i], tmp);
    fmpq_mat_mul(tmp, tmp, src);
  }

  fmpq_mat_init_set(pows[k], tmp);

  for (int i = 1; i <= m; i++) {
    fmpq_mat_mul(tmp, tmp, tmp);
    fmpq_mat_init_set(pows_exp[i], tmp);
  }

  fmpq_poly_t f_fmpq;
  fmpq_poly_init(f_fmpq);
  fmpq_poly_set_fmpz_poly(f_fmpq, f);

  _fmpq_poly_apply_fmpq_mat_ps_recur(dst, src, f_fmpq, pows, k, pows_exp, m);

  // fmpq_mat_print(dst);
  // printf("\n");

  for (int i = 0; i <= k + m; i++) {
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
  // else {
  //   fmpz_poly_apply_fmpq_mat_ps_recur(dst, src, f);
  //   return;
  // }
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



