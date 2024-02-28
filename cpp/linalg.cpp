#include "linalg.h"
#include "manin_basis.h"
#include "manin_element.h"

#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>

#include <vector>
#include <cassert>

std::vector<ManinElement> map_kernel(std::vector<ManinElement> B, std::function<ManinElement(ManinBasisElement)> f, int64_t M) {
  // If the input space is trivial, the result is also trivial.
  if (B.size() == 0) {
    return std::vector<ManinElement>();
  }

  // Otherwise, we find the level from the first element.
  int64_t N = B[0].N;
  std::vector<ManinBasisElement> N_basis = manin_basis(N);
  std::vector<ManinBasisElement> M_basis = manin_basis(M);

  // Construct matrix of B
  // TODO: just make a function that converts between sparse and dense reps of
  // a ManinElement basis
  fmpq_mat_t B_matrix;
  fmpq_mat_init(B_matrix, N_basis.size(), B.size());
  fmpq_mat_zero(B_matrix);

  fmpz_mat_t B_matrix_z;
  fmpz_mat_init(B_matrix_z, N_basis.size(), B.size());
  fmpz_mat_zero(B_matrix_z);

  for (int col = 0; col < B.size(); col++) {
    ManinElement b = B[col];
    // b.print_with_generators();
    // printf("\n");
    for (MBEWC component : b.components) {
      int row = component.basis_index;
      // printf("%d\n", row);
      assert(row < N_basis.size());
      fmpq_set(fmpq_mat_entry(B_matrix, row, col), &component.coeff);
    }
  }

  // XXX: this feels a bit wasteful, also is colwise right??
  fmpq_mat_get_fmpz_mat_colwise(B_matrix_z, NULL, B_matrix);
  fmpq_mat_clear(B_matrix);

  printf("basis matrix:\n");
  fmpz_mat_print_pretty(B_matrix_z);
  printf("\n");

  // Construct matrix of the map
  fmpq_mat_t map_matrix;
  fmpq_mat_init(map_matrix, M_basis.size(), B.size());
  fmpq_mat_zero(map_matrix);

  fmpz_mat_t map_matrix_z;
  fmpz_mat_init(map_matrix_z, M_basis.size(), B.size());
  fmpz_mat_zero(map_matrix_z);

  for (int col = 0; col < B.size(); col++) {
    // [ ]: maybe inline map() and cache f(mbe)?
    ManinElement fb = B[col].map(f, M);
    // fb.print_with_generators();
    // printf("\n");
    for (MBEWC component : fb.components) {
      int row = component.basis_index;
      // printf("%d\n", row);
      assert(row < M_basis.size());
      fmpq_set(fmpq_mat_entry(map_matrix, row, col), &component.coeff);
    }
  }

  // XXX: this feels a bit wasteful
  fmpq_mat_get_fmpz_mat_colwise(map_matrix_z, NULL, map_matrix);
  fmpq_mat_clear(map_matrix);

  // printf("map matrix:\n");
  // fmpz_mat_print_pretty(map_matrix_z);
  // printf("\n");

  fmpz_mat_t map_kernel, map_kernel_window;
  fmpz_mat_init(map_kernel, B.size(), B.size());
  int64_t rank = fmpz_mat_nullspace(map_kernel, map_matrix_z);
  fmpz_mat_window_init(map_kernel_window, map_kernel, 0, 0, B.size(), rank);

  fmpz_mat_clear(map_matrix_z);

  // printf("map kernel window:\n");
  // fmpz_mat_print_pretty(map_kernel_window);
  // printf("\n");

  fmpz_mat_t map_kernel_in_orig_basis;
  fmpz_t den;
  fmpz_init(den);
  fmpz_mat_init(map_kernel_in_orig_basis, N_basis.size(), rank);
  fmpz_mat_mul(map_kernel_in_orig_basis, B_matrix_z, map_kernel_window);
  fmpz_mat_content(den, map_kernel_in_orig_basis);
  if (!fmpz_is_zero(den)) {
    fmpz_mat_scalar_divexact_fmpz(map_kernel_in_orig_basis, map_kernel_in_orig_basis, den);
  }

  fmpz_mat_window_clear(map_kernel_window);
  fmpz_mat_clear(map_kernel);

  // printf("map kernel in orig basis:\n");
  // fmpz_mat_print_pretty(map_kernel_in_orig_basis);
  // printf("\n");


  // Convert each column of the kernel to a ManinElement
  std::vector<ManinElement> output;

  for (int col = 0; col < rank; col++) {
    std::vector<MBEWC> components;
    for (int row = 0; row < N_basis.size(); row++) {
      if (!(fmpz_is_zero(fmpz_mat_entry(map_kernel_in_orig_basis, row, col)))) {
        fmpq_t coeff;
        fmpq_init(coeff);
        fmpq_set_fmpz(coeff, fmpz_mat_entry(map_kernel_in_orig_basis, row, col));
        components.push_back({.basis_index = row, .coeff = *coeff});
      }
    }
    ManinElement element = {.N = N, .components = components};
    element.sort();
    output.push_back(element);
  }

  fmpz_mat_clear(map_kernel_in_orig_basis);

  return output;
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


void fmpz_poly_apply_fmpq_mat(fmpq_mat_t dst, const fmpq_mat_t src, const fmpz_poly_t f) {
  fmpz_t coeff;
  fmpz_init(coeff);
  int degree = fmpz_poly_degree(f);

  printf("[debug] fmpz_poly_apply_fmpq_mat called with "); // arguments\n");
  // printf("src: \n");
  // fmpq_mat_print(src);
  printf("f: ");
  fmpz_poly_print_pretty(f, "x");
  printf("\n");

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

DecomposeResult DecomposeResult::empty() {
  return {
    .done = std::vector<std::vector<ManinElement>>(),
    .remaining = std::vector<std::vector<ManinElement>>()
  };
}

// TODO: should be faster to decompose using a matrix of the action on the newform subspace

DecomposeResult decompose(std::vector<ManinElement> B, std::function<ManinElement(ManinBasisElement)> f) {
  // If the input space is trivial, the result is also trivial.
  if (B.size() == 0) {
    return DecomposeResult::empty();
  }

  // Otherwise, we find the level from the first element.
  int64_t N = B[0].N;
  std::vector<ManinBasisElement> N_basis = manin_basis(N);

  // Construct matrix of B
  // TODO: just make a function that converts between sparse and dense reps of
  // a ManinElement basis
  fmpq_mat_t B_matrix;
  fmpq_mat_init(B_matrix, N_basis.size(), B.size());
  fmpq_mat_zero(B_matrix);

  fmpz_mat_t B_matrix_z;
  fmpz_mat_init(B_matrix_z, N_basis.size(), B.size());
  fmpz_mat_zero(B_matrix_z);

  for (int col = 0; col < B.size(); col++) {
    ManinElement b = B[col];
    // b.print_with_generators();
    // printf("\n");
    for (MBEWC component : b.components) {
      int row = component.basis_index;
      // printf("%d\n", row);
      assert(row < N_basis.size());
      fmpq_set(fmpq_mat_entry(B_matrix, row, col), &component.coeff);
    }
  }

  // XXX: this feels a bit wasteful
  fmpq_mat_get_fmpz_mat_colwise(B_matrix_z, NULL, B_matrix);
  fmpq_mat_clear(B_matrix);

  // Finds pivot rows of the matrix B.
  std::vector<int> pivots;
  fmpz* pivot_coeffs = _fmpz_vec_init(B.size());

  printf("[debug] B_matrix_z:\n");
  fmpz_mat_print_pretty(B_matrix_z);
  printf("\n");

  int current_col = 0;
  for (int row = 0; row < N_basis.size(); row++) {
    bool is_pivot = true;
    for (int col = 0; col < B.size(); col++) {
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
  }

  for (int i = 0; i < B.size(); i++) {
    printf("%d ", pivots[i]);
    fmpz_print(pivot_coeffs + i);
    printf("\n");
  }

  // XXX: this section is here for debugging purposes
  // Construct matrix of the map
  int debug = 1;
  if (debug) {
    fmpq_mat_t map_matrix;
    fmpq_mat_init(map_matrix, N_basis.size(), B.size());
    fmpq_mat_zero(map_matrix);

    for (int col = 0; col < B.size(); col++) {
      // [ ]: maybe inline map() and cache f(mbe)?
      ManinElement fb = B[col].map(f, N);
      fb.print_with_generators();
      printf("\n");
      for (MBEWC component : fb.components) {
        int row = component.basis_index;
        // printf("%d\n", row);
        assert(row < N_basis.size());
        fmpq_set(fmpq_mat_entry(map_matrix, row, col), &component.coeff);
      }
    }

    printf("[debug] map_matrix:\n");
    fmpq_mat_print(map_matrix);
    printf("\n");

    fmpq_mat_clear(map_matrix);
  }

  // Construct matrix of the linear map f acting on B
  fmpq_mat_t f_matrix;
  fmpq_mat_init(f_matrix, B.size(), B.size());

  for (int col = 0; col < B.size(); col++) {
    ManinElement fb = B[col].map(f, N);

    // printf("rows: ");

    auto it = fb.components.begin();
    int pivot_index = 0;

    while (true) {
      if (it == fb.components.end()) break;
      if (pivot_index == pivots.size()) break;

      int row = it->basis_index;

      if (row > pivots[pivot_index]) {
        pivot_index++;
      } else if (row == pivots[pivot_index]) {
        fmpq_t coeff;
        fmpq_set(coeff, &it->coeff);
        fmpq_div_fmpz(coeff, coeff, (pivot_coeffs + pivot_index));
        fmpq_set(fmpq_mat_entry(f_matrix, pivot_index, col), coeff);
        fmpq_clear(coeff);
        it++;
        pivot_index++;
      } else if (row < pivots[pivot_index]) {
        it++;
      }
    }

    // for (MBEWC component : fb.components) {
    //   // printf("%d ", component.basis_index);
    //   if (component.basis_index == pivots[pivot_index]) {
    //     fmpq_t coeff;
    //     fmpq_set(coeff, &component.coeff);
    //     fmpq_div_fmpz(coeff, coeff, (pivot_coeffs + pivot_index));
    //     fmpq_set(fmpq_mat_entry(f_matrix, pivot_index, col), coeff);
    //     fmpq_clear(coeff);
    //   }
    //   else {
    //     while (pivot_index < pivots.size() && component.basis_index > pivots[pivot_index]) pivot_index++;
    //     if (pivot_index == pivots.size()) break;
    //   }
    // }
    // printf("\n");
  }

  printf("[debug] f_matrix:\n");
  fmpq_mat_print(f_matrix);
  printf("\n");

  fmpq_poly_t min_poly;
  fmpz_poly_t min_poly_z;
  fmpq_poly_init(min_poly);
  fmpz_poly_init(min_poly_z);

  fmpq_mat_minpoly(min_poly, f_matrix);
  // fmpq_mat_clear(f_matrix);

  fmpq_poly_get_numerator(min_poly_z, min_poly);
  fmpq_poly_clear(min_poly);

  fmpz_poly_factor_t min_poly_factored;
  fmpz_poly_factor_init(min_poly_factored);
  fmpz_poly_factor(min_poly_factored, min_poly_z);
  fmpz_poly_clear(min_poly_z);

  int num_factors = min_poly_factored->num;
  if (num_factors == 0) {
    fmpz_poly_factor_clear(min_poly_factored);
    return DecomposeResult::empty();
  }

  std::vector<std::vector<ManinElement>> done;
  std::vector<std::vector<ManinElement>> remaining;

  fmpq_mat_t poly_on_f_matrix;
  fmpz_mat_t poly_on_f_matrix_z, poly_mat_kernel;
  fmpz_t den;
  fmpq_mat_init(poly_on_f_matrix, B.size(), B.size());
  fmpz_mat_init(poly_on_f_matrix_z, B.size(), B.size());
  fmpz_mat_init(poly_mat_kernel, B.size(), B.size());
  fmpz_init(den);

  for (int i = 0; i < num_factors; i++) {
    fmpz_mat_t poly_mat_kernel_window, poly_mat_kernel_in_orig_basis;
    fmpz_poly_struct *factor = min_poly_factored->p + i;
    fmpz_poly_apply_fmpq_mat(poly_on_f_matrix, f_matrix, factor);
    int degree = fmpz_poly_degree(factor);
    fmpq_mat_get_fmpz_mat_colwise(poly_on_f_matrix_z, NULL, poly_on_f_matrix);
    int rank = fmpz_mat_nullspace(poly_mat_kernel, poly_on_f_matrix_z);
    fmpz_mat_window_init(poly_mat_kernel_window, poly_mat_kernel, 0, 0, B.size(), rank);

    fmpz_mat_init(poly_mat_kernel_in_orig_basis, N_basis.size(), rank);
    fmpz_mat_mul(poly_mat_kernel_in_orig_basis, B_matrix_z, poly_mat_kernel_window);
    fmpz_mat_window_clear(poly_mat_kernel_window);

    fmpz_mat_content(den, poly_mat_kernel_in_orig_basis);
    if (!fmpz_is_zero(den)) {
      fmpz_mat_scalar_divexact_fmpz(poly_mat_kernel_in_orig_basis, poly_mat_kernel_in_orig_basis, den);
    }

    std::vector<ManinElement> output;
    for (int col = 0; col < rank; col++) {
      std::vector<MBEWC> components;
      for (int row = 0; row < N_basis.size(); row++) {
        if (!(fmpz_is_zero(fmpz_mat_entry(poly_mat_kernel_in_orig_basis, row, col)))) {
          fmpq_t coeff;
          fmpq_init(coeff);
          fmpq_set_fmpz(coeff, fmpz_mat_entry(poly_mat_kernel_in_orig_basis, row, col));
          components.push_back({.basis_index = row, .coeff = *coeff});
        }
      }
      ManinElement element = {.N = N, .components = components};
      element.sort();
      output.push_back(element);
    }

    fmpz_mat_clear(poly_mat_kernel_in_orig_basis);

    if (degree == rank) {
      done.push_back(output);
    } else {
      remaining.push_back(output);
    }
  }

  fmpq_mat_clear(f_matrix);
  fmpq_mat_clear(poly_on_f_matrix);
  fmpz_mat_clear(poly_on_f_matrix_z);
  fmpz_mat_clear(poly_mat_kernel);
  fmpz_clear(den);

  return {.done = done, .remaining = remaining};
}