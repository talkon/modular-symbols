#include "manin_basis.h"
#include "subspace.h"
#include "newform_subspaces.h"
#include "fmpz_mat_helpers.h"
#include "debug_utils.h"

#include <flint/ulong_extras.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>

#include <cassert>

int Subspace::dimension() const {
  return basis.mat->c;
}

void Subspace::print(int index) const {
  // Label
  const char alphabet[] = "abcdefghijklmnopqrstuvwxyz";
  char newform_label[11];
  newform_label[10] = '\0';
  int p = 9;
  if (index == 0) {
    newform_label[p] = alphabet[0];
    p--;
  }
  for (; index > 0 && p >= 0; p--) {
    newform_label[p] = alphabet[index % 26];
    index /= 26;
  }
  if (p < 0) assert(false);
  printf("%d.2.a.%s:", level, newform_label + p + 1);

  // Dimension
  int dim = dimension();
  printf("%d:", dim);

  // Atkin-Lehner signs
  n_factor_t factors;
  n_factor_init(&factors);
  n_factor(&factors, level, 1);

  for (int i = 0; i < factors.num; i++) {
    int64_t p = factors.p[i];
    if (std::find(atkin_lehner_pos.begin(), atkin_lehner_pos.end(), p) != atkin_lehner_pos.end()) {
      printf("+%lld", p);
    } else {
      printf("-%lld", p);
    }
    if (i < factors.num - 1) printf(",");
  }

  // Trace form
  printf(":");
  for (int i = 1; i <= trace_depth; i++) {
    printf("%lld,", trace_form.at(i));
  }

  // Hecke field polynomial
  printf(":");
  if (hecke_field_poly.has_value() && dim <= 20) {
    auto& poly = hecke_field_poly.value().poly;
    for (int i = dim; i >= 0; i--) {
      fmpz_print(fmpz_poly_get_coeff_ptr(poly, i));
      printf(",");
    }
  } else if (dim == 1) {
    printf("1,0,");
  }

  // Hecke minimal polynomials
  printf(":");
  for (auto const& [p, mp] : hecke_min_polys) {
    printf("(");
    auto& poly = mp.poly;
    int deg = fmpz_poly_degree(poly);
    printf("%lld;%lld;", p, deg);
    for (int i = deg; i >= 0; i--) {
      fmpz_print(fmpz_poly_get_coeff_ptr(poly, i));
      printf(",");
    }
    printf(")");
  }


}

void Subspace::set_first_trace() {
  FmpqMatrix unused;
  next_trace(1, unused, -1);
}

int Subspace::next_trace(int next_depth, FmpqMatrix& hecke_mat, int max_trace_depth) {
  int n = trace_depth + 1;
  assert(n == next_depth);
  int dim = dimension();
  assert(dim > 0);

  n_factor_t factors;
  n_factor_init(&factors);
  n_factor(&factors, n, 1);

  int g = n_gcd(n, level);
  int is_prime = n_is_prime(n);

  if (g == 1) {
    fmpq_mat_t f_matrix;
    fmpq_mat_init(f_matrix, dim, dim);

    if (n == 1) {
      fmpq_mat_one(f_matrix);
    }
    else if (is_prime) {
      int p = n;
      // TODO: There's a lot of copied code here -- might be possible to refactor
      std::vector<ManinBasisElement> N_basis = manin_basis(level);

      // Finds pivot rows of the matrix B.
      std::vector<int> pivots;
      fmpz* pivot_coeffs = _fmpz_vec_init(dim);

      int current_col = 0;
      for (int row = 0; row < N_basis.size(); row++) {
        bool is_pivot = true;
        for (int col = 0; col < dim; col++) {
          if (col != current_col && !fmpz_is_zero(fmpz_mat_entry(basis.mat, row, col))) {
            is_pivot = false;
            break;
          }
        }

        if (!is_pivot) continue;

        if (fmpz_is_zero(fmpz_mat_entry(basis.mat, row, current_col))) continue;

        pivots.push_back(row);
        fmpz_set(pivot_coeffs + current_col, fmpz_mat_entry(basis.mat, row, current_col));
        current_col++;

        if (current_col == dim) break;
      }

      // Construct matrix of the linear map f acting on B
      fmpq_mat_t pivot_rows;
      fmpq_mat_init(pivot_rows, dim, N_basis.size());

      for (int row = 0; row < dim; row++) {
        for (int col = 0; col < N_basis.size(); col++) {
          fmpq_div_fmpz(
            fmpq_mat_entry(pivot_rows, row, col),
            fmpq_mat_entry(hecke_mat.mat, pivots[row], col),
            (pivot_coeffs + row)
          );
        }
      }

      DEBUG_INFO(6,
        {
          printf("pivot_rows: ");
          fmpq_mat_print_dimensions(pivot_rows);
          printf("\n");
        }
      )

      fmpq_mat_mul_fmpz_mat(f_matrix, pivot_rows, basis.mat);
      fmpq_mat_clear(pivot_rows);
      flint_cleanup();

      DEBUG_INFO(6,
        {
          printf("f_matrix: ");
          fmpq_mat_print_dimensions(f_matrix);
          printf("\n");
        }
      )

      _fmpz_vec_clear(pivot_coeffs, dim);
    }
    else if (factors.num == 1) {
      int64_t p = factors.p[0];
      int exp = factors.exp[0];
      int n_over_p = n_pow(p, exp - 1);
      fmpq_mat_mul(f_matrix, hecke_matrices.at(n_over_p).mat, hecke_matrices.at(p).mat);
      if (level % p != 0) {
        fmpq_mat_t temp;
        fmpz_t P;
        fmpq_mat_init(temp, dim, dim);
        fmpz_init_set_si(P, p);

        fmpq_mat_scalar_mul_fmpz(temp, hecke_matrices.at(n_over_p / p).mat, P);
        fmpq_mat_sub(f_matrix, f_matrix, temp);

        fmpq_mat_clear(temp);
        fmpz_clear(P);
        flint_cleanup();
      }
    }
    else {
      int64_t q = n_pow(factors.p[0], factors.exp[0]);
      fmpq_mat_set(f_matrix, hecke_matrices.at(q).mat);
      for (int i = 1; i < factors.num; i++) {
        q = n_pow(factors.p[i], factors.exp[i]);
        fmpq_mat_mul(f_matrix, f_matrix, hecke_matrices.at(q).mat);
      }
    }
    fmpq_t trace;
    fmpq_init(trace);
    fmpq_mat_trace(trace, f_matrix);
    assert(fmpz_is_one(fmpq_denref(trace)));

    int trace_int = fmpz_get_si(fmpq_numref(trace));
    trace_form.insert(std::make_pair(n, trace_int));
    fmpq_clear(trace);

    if (is_prime && dimension() <= 20) {
      fmpq_poly_t min_poly;
      fmpz_poly_t min_poly_z;
      fmpq_poly_init(min_poly);
      fmpz_poly_init(min_poly_z);
      if (fmpq_mat_is_zero(f_matrix)) {
        // Manually set the minimal polynomial to x.
        fmpz_poly_set_coeff_si(min_poly_z, 1, 1);
      } else {
        fmpq_mat_minpoly(min_poly, f_matrix);
        fmpq_poly_get_numerator(min_poly_z, min_poly);
      }

      fmpq_poly_clear(min_poly);

      FmpzPoly poly;
      poly.set_move(min_poly_z);
      hecke_min_polys.insert(std::make_pair(n, poly));
    }

    if (factors.num <= 1 && (max_trace_depth == -1 || 2 * n <= max_trace_depth)) {
      FmpqMatrix hecke_matrix;
      hecke_matrix.set_move(f_matrix);
      hecke_matrices.insert(std::make_pair(n, hecke_matrix));
    } else {
      fmpq_mat_clear(f_matrix);
      flint_cleanup();
    }

    trace_depth++;
  }
  else {
    // For primes p | N:
    // - if p^2 | N, then the action of the Hecke operator T_p is zero, and
    // - if p^2 \nmid N, then the action of the Hecke operator T_p is -W_p, where W_p is the Atkin-Lehner involution.
    // See Prop 13.3.4 in (Cohen, Stromberg)
    int64_t x = 1;
    int sign = 1;

    for (int i = 0; i < factors.num; i++) {
      int64_t p = factors.p[i];
      int exp = factors.exp[i];

      if (level % p == 0) {
        if (level % (p * p) == 0) sign = 0;
        else if (
          exp % 2 == 1
          && std::find(atkin_lehner_pos.begin(), atkin_lehner_pos.end(), p) != atkin_lehner_pos.end()
        ) {
          sign = -sign;
        }
      }
      else {
        x *= n_pow(p, exp);
      }
    }

    int trace_int = trace_form.at(x) * sign;
    trace_form.insert(std::make_pair(n, trace_int));
    trace_depth++;

    if (is_prime && dimension() < 20) {
      if (level % (n * n) == 0) {
        sign = 0;
      } else if (
        std::find(atkin_lehner_pos.begin(), atkin_lehner_pos.end(), n) != atkin_lehner_pos.end()
      ) {
        sign = -1;
      } else {
        sign = 1;
      }
      fmpz_poly_t min_poly;
      fmpz_poly_init(min_poly);
      fmpz_poly_set_coeff_si(min_poly, 1, 1);
      fmpz_poly_set_coeff_si(min_poly, 0, -sign);

      FmpzPoly poly;
      poly.set_move(min_poly);
      hecke_min_polys.insert(std::make_pair(n, poly));
    }

  }

  return trace_depth;
}

void Subspace::clear_hecke_matrices() {
  hecke_matrices.clear();
}