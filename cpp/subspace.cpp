#include "manin_basis.h"
#include "subspace.h"
#include "newform_subspaces.h"
#include "fmpz_mat_helpers.h"
#include "debug_utils.h"

#include <flint/ulong_extras.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>

#include <cassert>

int Subspace::dimension() const {
  return basis.size();
}

void Subspace::print() const {
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
  int64_t level = basis[0].N;

  n_factor_t factors;
  n_factor_init(&factors);
  n_factor(&factors, n, 1);

  int g = n_gcd(n, level);

  if (g == 1) {
    fmpq_mat_t f_matrix;
    fmpq_mat_init(f_matrix, dim, dim);

    if (n == 1) {
      fmpq_mat_one(f_matrix);
    }
    else if (n_is_prime(n)) {
      int p = n;
      // TODO: There's a lot of copied code here -- might be possible to refactor
      std::vector<ManinBasisElement> N_basis = manin_basis(level);

      // Construct matrix of B
      fmpq_mat_t B_matrix;
      fmpq_mat_init(B_matrix, N_basis.size(), basis.size());
      fmpq_mat_zero(B_matrix);

      fmpz_mat_t B_matrix_z;
      fmpz_mat_init(B_matrix_z, N_basis.size(), basis.size());
      fmpz_mat_zero(B_matrix_z);

      for (int col = 0; col < basis.size(); col++) {
        ManinElement b = basis[col];
        for (MBEWC component : b.components) {
          int row = component.basis_index;
          assert(row < N_basis.size());
          fmpq_set(fmpq_mat_entry(B_matrix, row, col), component.coeff);
        }
      }

      fmpq_mat_get_fmpz_mat_colwise(B_matrix_z, NULL, B_matrix);
      fmpq_mat_clear(B_matrix);

      DEBUG_INFO(6,
        {
          printf("B_matrix_z: ");
          fmpz_mat_print_dimensions(B_matrix_z);
          printf("\n");
        }
      )

      // Finds pivot rows of the matrix B.
      std::vector<int> pivots;
      fmpz* pivot_coeffs = _fmpz_vec_init(basis.size());

      int current_col = 0;
      for (int row = 0; row < N_basis.size(); row++) {
        bool is_pivot = true;
        for (int col = 0; col < basis.size(); col++) {
          if (col != current_col && !fmpz_is_zero(fmpz_mat_entry(B_matrix_z, row, col))) {
            is_pivot = false;
            break;
          }
        }

        if (!is_pivot) continue;

        if (fmpz_is_zero(fmpz_mat_entry(B_matrix_z, row, current_col))) continue;

        pivots.push_back(row);
        fmpz_set(pivot_coeffs + current_col, fmpz_mat_entry(B_matrix_z, row, current_col));
        current_col++;

        if (current_col == basis.size()) break;
      }

      // Construct matrix of the linear map f acting on B
      fmpq_mat_t pivot_rows;
      fmpq_mat_init(pivot_rows, basis.size(), N_basis.size());

      for (int row = 0; row < basis.size(); row++) {
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

      fmpq_mat_mul_fmpz_mat(f_matrix, pivot_rows, B_matrix_z);
      fmpq_mat_clear(pivot_rows);
      fmpz_mat_clear(B_matrix_z);

      DEBUG_INFO(6,
        {
          printf("f_matrix: ");
          fmpq_mat_print_dimensions(f_matrix);
          printf("\n");
        }
      )

      _fmpz_vec_clear(pivot_coeffs, basis.size());
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

    if (factors.num <= 1 && (max_trace_depth == -1 || 2 * n < max_trace_depth)) {
      FmpqMatrix hecke_matrix;
      hecke_matrix.set_move(f_matrix);
      hecke_matrices.insert(std::make_pair(n, hecke_matrix));
    } else {
      fmpq_mat_clear(f_matrix);
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
  }

  return trace_depth;
}

void Subspace::clear_hecke_matrices() {
  hecke_matrices.clear();
}