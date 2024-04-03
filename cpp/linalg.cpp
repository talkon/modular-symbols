#include "linalg.h"
#include "manin_basis.h"
#include "manin_element.h"
#include "debug_utils.h"
#include "fmpz_mat_helpers.h"

#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>

#include <vector>
#include <cassert>
#include <stdexcept>

// TODO: consider requiring f to return a reference?
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
      fmpq_set(fmpq_mat_entry(B_matrix, row, col), component.coeff);
    }
  }

  // XXX: this feels a bit wasteful, also is colwise right??
  fmpq_mat_get_fmpz_mat_colwise(B_matrix_z, NULL, B_matrix);
  fmpq_mat_clear(B_matrix);

  DEBUG_INFO(4,
    {
      printf("B_matrix_z: ");
      fmpz_mat_print_dimensions(B_matrix_z);
      printf("\n");
    }
  )

  // Construct matrix of the map
  fmpq_mat_t map_matrix;
  fmpq_mat_init(map_matrix, M_basis.size(), B.size());
  fmpq_mat_zero(map_matrix);

  fmpz_mat_t map_matrix_z;
  fmpz_mat_init(map_matrix_z, M_basis.size(), B.size());
  fmpz_mat_zero(map_matrix_z);

  bool use_map_of_basis = true;

  if (use_map_of_basis) {
    fmpq_mat_t map_of_basis;
    fmpq_mat_init(map_of_basis, M_basis.size(), N_basis.size());
    for (int col = 0; col < N_basis.size(); col++) {
      ManinElement fb = f(N_basis[col]);

      for (MBEWC component : fb.components) {
        int row = component.basis_index;
        // printf("%d\n", row);
        assert(row < M_basis.size());
        fmpq_set(fmpq_mat_entry(map_of_basis, row, col), component.coeff);
      }
    }

    fmpq_mat_mul_fmpz_mat(map_matrix, map_of_basis, B_matrix_z);
    fmpq_mat_clear(map_of_basis);

  } else {

    for (int col = 0; col < B.size(); col++) {
      // [ ]: maybe inline map()?
      ManinElement fb = B[col].map(f, M);
      for (MBEWC component : fb.components) {
        int row = component.basis_index;
        // printf("%d\n", row);
        assert(row < M_basis.size());
        fmpq_set(fmpq_mat_entry(map_matrix, row, col), component.coeff);
      }
    }
  }

  // XXX: this feels a bit wasteful
  fmpq_mat_get_fmpz_mat_rowwise(map_matrix_z, NULL, map_matrix);
  fmpq_mat_clear(map_matrix);

  // fflush(stdout);

  fmpz_mat_t map_kernel, map_kernel_window;
  // printf("aa %zu\n", B.size());
  fmpz_mat_init(map_kernel, B.size(), B.size());
  // printf("bb %zu\n", B.size());
  DEBUG_INFO(4,
    {
      printf("map_matrix_z: ");
      fmpz_mat_print_dimensions(map_matrix_z);
      printf("\n");
    }
  )

  // TODO: consider forcing this to use rref_mul instead of rref_fflu
  int64_t rank = fmpz_mat_nullspace_mul(map_kernel, map_matrix_z);
  // printf("cc %lld\n", rank);

  if (false) {
    fmpz_t d;
    fmpz_init(d);

    fmpz_mat_t tmp;
    fmpz_mat_init(tmp, M_basis.size(), B.size());

    DEBUG_INFO_PRINT(4, "computing tmp\n");
    fmpz_mat_rref_fflu(tmp, d, map_matrix_z);

    DEBUG_INFO(4,
      {
        printf("tmp [fflu]: ");
        fmpz_mat_print_dimensions(tmp);
        printf("\n");
      }
    )

    fmpz_mat_div_colwise_gcd(tmp);

    DEBUG_INFO(4,
      {
        printf("tmp: ");
        fmpz_mat_print_dimensions(tmp);
        printf("\n");
      }
    )

    fmpz_mat_clear(tmp);

    fmpz_mat_t tmp2;
    fmpz_mat_init(tmp2, M_basis.size(), B.size());

    DEBUG_INFO_PRINT(4, "computing tmp2\n");
    fmpz_mat_rref_mul(tmp2, d, map_matrix_z);

    DEBUG_INFO(4,
      {
        printf("tmp2 [mul]: ");
        fmpz_mat_print_dimensions(tmp2);
        printf("\n");
      }
    )

    fmpz_mat_div_colwise_gcd(tmp2);

    DEBUG_INFO(4,
      {
        printf("tmp2: ");
        fmpz_mat_print_dimensions(tmp2);
        printf("\n");
      }
    )

    fmpz_mat_clear(tmp2);
    fmpz_clear(d);
  }

  fmpz_mat_window_init(map_kernel_window, map_kernel, 0, 0, B.size(), rank);

  DEBUG_INFO(4,
    {
      printf("kernel: ");
      fmpz_mat_print_dimensions(map_kernel_window);
      printf("\n");
    }
  )

  fmpz_mat_div_colwise_gcd(map_kernel_window);

  fmpz_mat_clear(map_matrix_z);

  // printf("map kernel window:\n");
  // fmpz_mat_print_pretty(map_kernel_window);
  // printf("\n");

  fmpz_mat_t map_kernel_in_orig_basis;
  fmpz_mat_init(map_kernel_in_orig_basis, N_basis.size(), rank);

  DEBUG_INFO(4,
    {
      printf("kernel (cleared): ");
      fmpz_mat_print_dimensions(map_kernel_window);
      printf("\n");
    }
  )

  // XXX: This multiplication is still slow, maybe force B_matrix_z and map_kernel_window to be in rref form?
  fmpz_mat_mul(map_kernel_in_orig_basis, B_matrix_z, map_kernel_window);

  DEBUG_INFO(4,
    {
      printf("result: ");
      fmpz_mat_print_dimensions(map_kernel_in_orig_basis);
      printf("\n");
    }
  )

  // printf("dd %lld\n", rank);

  fmpz_mat_clear(B_matrix_z);

  fmpz_mat_div_colwise_gcd(map_kernel_in_orig_basis);

  fmpz_mat_window_clear(map_kernel_window);
  fmpz_mat_clear(map_kernel);

  DEBUG_INFO(4,
    {
      printf("result (cleared): ");
      fmpz_mat_print_dimensions(map_kernel_in_orig_basis);
      printf("\n");
    }
  )


  // printf("map kernel in orig basis:\n");
  // fmpz_mat_print_pretty(map_kernel_in_orig_basis);
  // printf("\n");

  // Convert each column of the kernel to a ManinElement
  std::vector<ManinElement> output;

  for (int col = 0; col < rank; col++) {
    // printf("col %d\n", col);
    std::vector<MBEWC> components;
    for (int row = 0; row < N_basis.size(); row++) {
      if (!(fmpz_is_zero(fmpz_mat_entry(map_kernel_in_orig_basis, row, col)))) {
        fmpq_t coeff;
        fmpq_init(coeff);
        fmpq_set_fmpz(coeff, fmpz_mat_entry(map_kernel_in_orig_basis, row, col));
        components.emplace_back(row, coeff);
        fmpq_clear(coeff);
      }
    }
    has_duplicate_keys(components);
    ManinElement element = ManinElement(N, components);
    element.mark_as_sorted_unchecked();
    output.push_back(element);
  }

  // printf("ff %lld\n", rank);

  // printf("output matrix:\n");
  // fmpz_mat_print_pretty(map_kernel_in_orig_basis);
  // printf("\n");

  fmpz_mat_clear(map_kernel_in_orig_basis);

  return output;
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
      fmpq_set(fmpq_mat_entry(B_matrix, row, col), component.coeff);
    }
  }

  // printf("[debug] B_matrix:\n");
  // fmpq_mat_print(B_matrix);
  // printf("\n");

  // XXX: this feels a bit wasteful
  fmpq_mat_get_fmpz_mat_colwise(B_matrix_z, NULL, B_matrix);
  fmpq_mat_clear(B_matrix);

  // Finds pivot rows of the matrix B.
  std::vector<int> pivots;
  fmpz* pivot_coeffs = _fmpz_vec_init(B.size());

  // printf("[debug] B_matrix_z:\n");
  // fmpz_mat_print_pretty(B_matrix_z);
  // printf("\n");

  int current_col = 0;
  for (int row = 0; row < N_basis.size(); row++) {
    // printf("%d %d\n", row, current_col);
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

    if (current_col == B.size()) break;
  }

  // for (int i = 0; i < B.size(); i++) {
  //   printf("%d ", pivots[i]);
  //   fmpz_print(pivot_coeffs + i);
  //   printf("\n");
  // }

  // Construct matrix of the linear map f acting on B
  fmpq_mat_t f_matrix;
  fmpq_mat_init(f_matrix, B.size(), B.size());

  bool use_map_of_basis = true;

  if (use_map_of_basis) {
    fmpq_mat_t map_of_basis;
    fmpq_mat_init(map_of_basis, B.size(), N_basis.size());
    for (int col = 0; col < N_basis.size(); col++) {
      ManinElement fb = f(N_basis[col]);

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
          fmpq_init(coeff);
          fmpq_set(coeff, it->coeff);
          fmpq_div_fmpz(coeff, coeff, (pivot_coeffs + pivot_index));
          fmpq_set(fmpq_mat_entry(map_of_basis, pivot_index, col), coeff);
          fmpq_clear(coeff);
          it++;
          pivot_index++;
        } else if (row < pivots[pivot_index]) {
          it++;
        }
      }
    }

    fmpq_mat_mul_fmpz_mat(f_matrix, map_of_basis, B_matrix_z);
    fmpq_mat_clear(map_of_basis);

  } else {
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
          fmpq_init(coeff);
          fmpq_set(coeff, it->coeff);
          fmpq_div_fmpz(coeff, coeff, (pivot_coeffs + pivot_index));
          fmpq_set(fmpq_mat_entry(f_matrix, pivot_index, col), coeff);
          fmpq_clear(coeff);
          it++;
          pivot_index++;
        } else if (row < pivots[pivot_index]) {
          it++;
        }
      }
    }
  }

  _fmpz_vec_clear(pivot_coeffs, B.size());

  fmpq_poly_t min_poly;
  fmpz_poly_t min_poly_z;
  fmpq_poly_init(min_poly);
  fmpz_poly_init(min_poly_z);

  // XXX: FLINT seems to think that the minimal polynomial of the zero matrix is 1, and not T.
  if (fmpq_mat_is_zero(f_matrix)) {
    // Manually set the minimal polynomial to x.
    fmpz_poly_set_coeff_si(min_poly_z, 1, 1);
  } else {
    fmpq_mat_minpoly(min_poly, f_matrix);
    fmpq_poly_get_numerator(min_poly_z, min_poly);
  }

  DEBUG_INFO_PRINT(3, " min_poly_z degree: %ld\n", fmpz_poly_degree(min_poly_z));

  DEBUG_INFO(4,
    {
      printf(" min_poly_z: ");
      fmpz_poly_print_pretty(min_poly_z, "T");
      printf("\n");
    }
  )

  fmpq_poly_clear(min_poly);

  fmpz_poly_factor_t min_poly_factored;
  fmpz_poly_factor_init(min_poly_factored);
  fmpz_poly_factor(min_poly_factored, min_poly_z);

  std::vector<std::vector<ManinElement>> done;
  std::vector<std::vector<ManinElement>> remaining;

  int num_factors = min_poly_factored->num;
  if (num_factors == 0) {
    // this should actually be impossible?
    assert(false);
  } else if (num_factors == 1 && fmpz_poly_degree(min_poly_z) == B.size()) {
    done.push_back(B);
    DEBUG_INFO(3,
      {
        printf(" minimal polynomial irreducible and degree equal to space dimension %zu\n", B.size());
      }
    )
  } else {
    fmpq_mat_t poly_on_f_matrix;
    fmpz_mat_t poly_on_f_matrix_z, poly_mat_kernel;
    fmpq_mat_init(poly_on_f_matrix, B.size(), B.size());
    fmpz_mat_init(poly_on_f_matrix_z, B.size(), B.size());
    fmpz_mat_init(poly_mat_kernel, B.size(), B.size());

    for (int i = 0; i < num_factors; i++) {
      fmpz_mat_t poly_mat_kernel_window, poly_mat_kernel_in_orig_basis;
      fmpz_poly_struct *factor = min_poly_factored->p + i;

      // TODO: when deg(factor) is large (>6?) it's probably faster to first compute the subspace (after modding out by small factors)
      // and verify that the minimal polynomial on that subspace has max degree.
      fmpz_poly_apply_fmpq_mat_ps(poly_on_f_matrix, f_matrix, factor);

      int degree = fmpz_poly_degree(factor);
      // NOTE: this needs to be rowwise!
      fmpq_mat_get_fmpz_mat_rowwise(poly_on_f_matrix_z, NULL, poly_on_f_matrix);
      int rank = fmpz_mat_nullspace_mul(poly_mat_kernel, poly_on_f_matrix_z);
      fmpz_mat_window_init(poly_mat_kernel_window, poly_mat_kernel, 0, 0, B.size(), rank);

      fmpz_mat_div_colwise_gcd(poly_mat_kernel_window);

      fmpz_mat_init(poly_mat_kernel_in_orig_basis, N_basis.size(), rank);
      fmpz_mat_mul(poly_mat_kernel_in_orig_basis, B_matrix_z, poly_mat_kernel_window);
      fmpz_mat_window_clear(poly_mat_kernel_window);

      fmpz_mat_div_colwise_gcd(poly_mat_kernel_in_orig_basis);

      std::vector<ManinElement> output;
      for (int col = 0; col < rank; col++) {
        std::vector<MBEWC> components;
        for (int row = 0; row < N_basis.size(); row++) {
          if (!(fmpz_is_zero(fmpz_mat_entry(poly_mat_kernel_in_orig_basis, row, col)))) {
            fmpq_t coeff;
            fmpq_init(coeff);
            fmpq_set_fmpz(coeff, fmpz_mat_entry(poly_mat_kernel_in_orig_basis, row, col));
            components.push_back(MBEWC(row, coeff));
            fmpq_clear(coeff);
          }
        }
        ManinElement element = ManinElement(N, components);
        element.mark_as_sorted_unchecked();
        output.push_back(element);
      }

      fmpz_mat_clear(poly_mat_kernel_in_orig_basis);

      DEBUG_INFO(3,
        {
          printf(" subspace dimension: %d, factor degree: %d\n", rank, degree);
        }
      )

      if (degree == rank) {
        done.push_back(output);
      } else {
        remaining.push_back(output);
      }
    }

    fmpq_mat_clear(poly_on_f_matrix);
    fmpz_mat_clear(poly_on_f_matrix_z);
    fmpz_mat_clear(poly_mat_kernel);
  }

  fmpz_poly_clear(min_poly_z);
  fmpq_mat_clear(f_matrix);
  fmpz_mat_clear(B_matrix_z);
  fmpz_poly_factor_clear(min_poly_factored);

  return {.done = done, .remaining = remaining};
}