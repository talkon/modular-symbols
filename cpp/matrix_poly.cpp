#include "matrix_poly.h"
#include "fmpz_mat_helpers.h"
#include "debug_utils.h"

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>
#include <flint/profiler.h>

#include <cassert>

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
void fmpz_poly_apply_fmpq_mat_ps(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f, const slong mem_threshold) {

  DEBUG_INFO(5,
    {
      printf("fmpz_poly_apply_fmpq_mat_ps called with f(T) = ");
      fmpz_poly_print_pretty(f, "T");
      printf("\n");
      printf("mem_threshold = %ld\n", mem_threshold);
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
  int m = 1;
  while (m < k) {

    slong mem_usage = -2;
    if (mem_threshold > 0) {
      meminfo_t meminfo;
      get_memory_usage(meminfo);
      mem_usage = meminfo->rss;
      DEBUG_INFO_PRINT(5, "fmpz_poly_apply_fmpq_mat_ps: m = %d, mem usage = %d\n", m, (int) meminfo->rss);
    }

    if (m <= 10 || mem_usage < mem_threshold) {
      fmpq_mat_init_set(pows[m], tmp);
      fmpq_mat_mul(tmp, tmp, src);
      m++;
    } else {
      DEBUG_INFO_PRINT(3, "fmpz_poly_apply_fmpq_mat_ps: stopped early at m = %d due to memory usage limit\n", m);
      break;
    }
  }

  fmpq_mat_init_set(pows[m], tmp);

  fmpq_mat_zero(dst);

  fmpz_t a;

  fmpz_init(a);
  for (int j = l / m; j >= 0; j--) {

    for (int i = 0; i < m; i++) {
      fmpz_poly_get_coeff_fmpz(a, f, j * m + i);
      if (!fmpz_is_zero(a)) {
        fmpq_mat_scalar_mul_fmpz(tmp, pows[i], a);
        fmpq_mat_add(dst, dst, tmp);
      }
    }

    if (mem_threshold > 0) {
      meminfo_t meminfo;
      get_memory_usage(meminfo);
      DEBUG_INFO_PRINT(5, "fmpz_poly_apply_fmpq_mat_ps: j = %d, mem usage = %d\n", j, (int) meminfo->rss);
    }

    if (j > 0) {
      fmpq_mat_mul(dst, dst, pows[m]);
    }
  }

  fmpz_clear(a);
  fmpq_mat_clear(tmp);
  for (int i = 0; i <= m; i++) {
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

  // fmpz_mat_print_pretty(dst_z);
  // printf("\nden_pow: ");
  // fmpz_print(den_pow);
  // printf("\n");

  fmpz_poly_clear(f_adj);
  fmpz_clear(den);
  fmpz_clear(den_pow);
  fmpz_mat_clear(src_z);
  fmpz_mat_clear(dst_z);
}

// Algorithm B, but does the whole computation mod a prime, then combining once at the end.
void fmpz_poly_apply_nmod_mat_ps(nmod_mat_t dst, const nmod_mat_t src, const fmpz_poly_t f) {

  ulong n = src->mod.n;

  DEBUG_INFO(5,
    {
      printf("fmpz_poly_apply_nmod_mat_ps called with modulus %lu\n", n);
    }
  )

  ulong l = fmpz_poly_degree(f);
  ulong k = n_sqrt(l);

  ulong d = nmod_mat_nrows(src);
  assert(d == nmod_mat_nrows(src));

  // Compute 1, T, T^2, .., T^k
  nmod_mat_t* pows = (nmod_mat_t*) flint_malloc((k + 1) * sizeof(nmod_mat_t));
  nmod_mat_t tmp;

  nmod_mat_init(pows[0], d, d, n);
  for (int r = 0; r < d; r++) {
    nmod_mat_set_entry(pows[0], r, r, 1);
  }

  nmod_mat_init_set(tmp, src);
  for (int i = 1; i < k; i++) {
    nmod_mat_init_set(pows[i], tmp);
    nmod_mat_mul(tmp, tmp, src);
  }

  nmod_mat_init_set(pows[k], tmp);

  nmod_mat_zero(dst);

  fmpz_t a;

  fmpz_init(a);
  for (int j = l / k; j >= 0; j--) {

    for (int i = 0; i < k; i++) {
      fmpz_poly_get_coeff_fmpz(a, f, j * k + i);
      if (!fmpz_is_zero(a)) {
        nmod_mat_scalar_mul_fmpz(tmp, pows[i], a);
        nmod_mat_add(dst, dst, tmp);
      }
    }

    if (j > 0) {
      nmod_mat_mul(dst, dst, pows[k]);
    }
  }

  fmpz_clear(a);
  nmod_mat_clear(tmp);
  for (int i = 0; i <= k; i++) {
    nmod_mat_clear(pows[i]);
  }
  flint_free(pows);
}

void fmpz_poly_apply_fmpq_mat_ps_nmod(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f) {

  DEBUG_INFO(5,
    {
      printf("fmpz_poly_apply_fmpq_mat_ps_nmod called with f(T) = ");
      fmpz_poly_print_pretty(f, "T");
      printf("\n");
    }
  )

  ulong d = fmpq_mat_nrows(src);
  assert(d == fmpq_mat_ncols(src));

  fmpz_mat_t src_z, dst_z;
  fmpz_t den, mod;
  fmpz_t min_elt, max_elt;

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

  fmpz_init_set_ui(mod, 1);
  fmpz_init_set_si(min_elt, 0);
  fmpz_init_set_si(max_elt, 0);

  ulong p = 1ULL << 61;
  int unchanged_iters = 0;
  const int iters_threshold = 1;
  bool done = false;

  for (int iter = 0; unchanged_iters < iters_threshold; iter++) {
    p = n_nextprime(p, 1);

    // fmpz_mat_print_pretty(dst_z);
    // printf("\n");

    nmod_mat_t mat_p;
    nmod_mat_init(mat_p, d, d, p);
    fmpz_mat_get_nmod_mat(mat_p, src_z);

    fmpz_poly_apply_nmod_mat_ps(mat_p, mat_p, f_adj);

    if (iter == 0) {
      fmpz_mat_set_nmod_mat(dst_z, mat_p);
    } else {
      fmpz_mat_CRT_ui(dst_z, dst_z, mod, mat_p, 1);
    }

    fmpz_t new_min_elt, new_max_elt;
    fmpz_init(new_min_elt);
    fmpz_init(new_max_elt);
    fmpz_mat_min_elt(new_min_elt, dst_z);
    fmpz_mat_max_elt(new_max_elt, dst_z);

    if (fmpz_equal(new_min_elt, min_elt) && fmpz_equal(new_max_elt, max_elt)) {
      unchanged_iters++;
    } else {
      unchanged_iters = 0;
      fmpz_set(min_elt, new_min_elt);
      fmpz_set(max_elt, new_max_elt);
    }

    fmpz_clear(new_min_elt);
    fmpz_clear(new_max_elt);

    fmpz_mul_ui(mod, mod, p);
    nmod_mat_clear(mat_p);
  }

  fmpq_mat_set_fmpz_mat_div_fmpz(dst, dst_z, den_pow);

  // fmpz_mat_print_pretty(dst_z);
  // printf("\nden_pow: ");
  // fmpz_print(den_pow);
  // printf("\n");

  fmpz_poly_clear(f_adj);
  fmpz_clear(mod);
  fmpz_clear(den);
  fmpz_clear(den_pow);
  fmpz_mat_clear(src_z);
  fmpz_mat_clear(dst_z);
}

// Algorithm C in Paterson-Stockmeyer: sqrt(2 * deg(P)) matrix multiplications, but more scalar preprocessing
// BUG: using this function seems to give wrong results in some cases, and is also just very slow!
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

// Main function
void fmpz_poly_apply_fmpq_mat(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f, const slong mem_threshold) {
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
  } else if (degree <= 3) {
    fmpz_poly_apply_fmpq_mat_horner(dst, src, f);
  } else {
    fmpz_poly_apply_fmpq_mat_ps(dst, src, f, mem_threshold);
    // fmpz_poly_apply_fmpq_mat_ps_nmod(dst, src, f);
    // fmpz_poly_apply_fmpq_mat_ps_clear_denom(dst, src, f);
  }
  // else {
  //   fmpz_poly_apply_fmpq_mat_ps_recur(dst, src, f);
  //   return;
  // }

  DEBUG_INFO(5,
    {
      printf("result = ");
      fmpq_mat_print(dst);
      printf("\n");
    }
  )
}