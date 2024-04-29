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

  // Convert each column of the kernel to a ManinElement
  std::vector<ManinElement> output;

  for (int col = 0; col < rank; col++) {
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

  fmpz_mat_clear(map_kernel_in_orig_basis);

  return output;
}

SplitResult SplitResult::empty() {
  return {
    .pos_space = std::vector<ManinElement>(),
    .neg_space = std::vector<ManinElement>()
  };
}

SplitResult split(std::vector<ManinElement> B, std::function<ManinElement (ManinBasisElement)> f) {
  // If the input space is trivial, the result is also trivial.
  if (B.size() == 0) {
    return SplitResult::empty();
  }

  // Otherwise, we find the level from the first element.
  int64_t N = B[0].N;
  std::vector<ManinBasisElement> N_basis = manin_basis(N);

  // Construct matrix of B
  // TODO: store B as a matrix and pass the matrix in instead.
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

  // XXX: this feels a bit wasteful
  fmpq_mat_get_fmpz_mat_colwise(B_matrix_z, NULL, B_matrix);
  fmpq_mat_clear(B_matrix);

  DEBUG_INFO(4,
    {
      printf("B_matrix_z: ");
      fmpz_mat_print_dimensions(B_matrix_z);
      printf("\n");
    }
  )

  // Finds pivot rows of the matrix B.
  std::vector<int> pivots;
  fmpz* pivot_coeffs = _fmpz_vec_init(B.size());

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

  // Construct matrix of the linear map f acting on B
  fmpq_mat_t f_matrix;
  fmpq_mat_init(f_matrix, B.size(), B.size());

  bool use_map_of_basis = true;

  if (use_map_of_basis) {
    // TODO: just pass in this matrix instead of computing it multiple times.
    fmpq_mat_t map_of_basis;
    fmpq_mat_init(map_of_basis, B.size(), N_basis.size());
    for (int col = 0; col < N_basis.size(); col++) {
      ManinElement fb = f(N_basis[col]);

      // printf("%zu out of %zu\n", fb.components.size(), N_basis.size());

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

    DEBUG_INFO(4,
      {
        printf("map_of_basis: ");
        fmpq_mat_print_dimensions(map_of_basis);
        printf("\n");
      }
    )

    fmpq_mat_mul_fmpz_mat(f_matrix, map_of_basis, B_matrix_z);
    fmpq_mat_clear(map_of_basis);

    DEBUG_INFO(4,
      {
        printf("f_matrix: ");
        fmpq_mat_print_dimensions(f_matrix);
        printf("\n");
      }
    )

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

  // Case 1: Space doesn't split: only +1 space
  if (fmpq_mat_is_one(f_matrix)) {
    DEBUG_INFO_PRINT(3, " space doesn't split: A-L action is +1\n");

    fmpq_mat_clear(f_matrix);
    fmpz_mat_clear(B_matrix_z);

    return {.pos_space = B, .neg_space = std::vector<ManinElement>()};
  }

  fmpq_mat_t neg_one;
  fmpq_mat_init(neg_one, B.size(), B.size());
  fmpq_mat_one(neg_one);
  fmpq_mat_neg(neg_one, neg_one);

  // Case 2: Space doesn't split: only -1 space
  if (fmpq_mat_equal(f_matrix, neg_one)) {
    DEBUG_INFO_PRINT(3, " space doesn't split: A-L action is -1\n");

    fmpq_mat_clear(f_matrix);
    fmpz_mat_clear(B_matrix_z);
    fmpq_mat_clear(neg_one);

    return {.pos_space = std::vector<ManinElement>(), .neg_space = B};
  }

  // Case 3: Space splits into +1 and -1 space
  fmpz_mat_t f_matrix_z, kernel, kernel_window, kernel_in_orig_basis;
  fmpz_mat_init(f_matrix_z, B.size(), B.size());
  fmpz_mat_init(kernel, B.size(), B.size());

  // Set f_matrix to T-1, i.e. positive space
  // XXX: should be faster to just add to the diagonal elements
  fmpq_mat_add(f_matrix, f_matrix, neg_one);
  fmpq_mat_get_fmpz_mat_rowwise(f_matrix_z, NULL, f_matrix);

  int rank = fmpz_mat_nullspace_mul(kernel, f_matrix_z);
  fmpz_mat_window_init(kernel_window, kernel, 0, 0, B.size(), rank);

  DEBUG_INFO(4,
    {
      printf("(+) kernel_window: ");
      fmpz_mat_print_dimensions(kernel_window);
      printf("\n");
    }
  )

  fmpz_mat_div_colwise_gcd(kernel_window);

  DEBUG_INFO(4,
    {
      printf("(+) kernel_window (cleared): ");
      fmpz_mat_print_dimensions(kernel_window);
      printf("\n");
    }
  )

  fmpz_mat_init(kernel_in_orig_basis, N_basis.size(), rank);
  fmpz_mat_mul(kernel_in_orig_basis, B_matrix_z, kernel_window);

  DEBUG_INFO(4,
    {
      printf("(+) kernel_in_orig_basis: ");
      fmpz_mat_print_dimensions(kernel_in_orig_basis);
      printf("\n");
    }
  )

  fmpz_mat_div_colwise_gcd(kernel_in_orig_basis);

  DEBUG_INFO(4,
    {
      printf("(+) kernel_in_orig_basis (cleared): ");
      fmpz_mat_print_dimensions(kernel_in_orig_basis);
      printf("\n");
    }
  )

  // TODO: we should just return the pos_space basis as a matrix instead
  // of this conversion to sparse rep.
  std::vector<ManinElement> pos_space;
  for (int col = 0; col < rank; col++) {
    std::vector<MBEWC> components;
    for (int row = 0; row < N_basis.size(); row++) {
      if (!(fmpz_is_zero(fmpz_mat_entry(kernel_in_orig_basis, row, col)))) {
        fmpq_t coeff;
        fmpq_init(coeff);
        fmpq_set_fmpz(coeff, fmpz_mat_entry(kernel_in_orig_basis, row, col));
        components.push_back(MBEWC(row, coeff));
        fmpq_clear(coeff);
      }
    }
    ManinElement element = ManinElement(N, components);
    element.mark_as_sorted_unchecked();
    pos_space.push_back(element);
  }

  DEBUG_INFO_PRINT(3, " pos_space dimension: %d\n", rank);

  // Set f_matrix from T-1 to T+1, i.e. positive space
  // XXX: should be faster to just add to the diagonal elements
  fmpq_mat_sub(f_matrix, f_matrix, neg_one);
  fmpq_mat_sub(f_matrix, f_matrix, neg_one);
  fmpq_mat_get_fmpz_mat_rowwise(f_matrix_z, NULL, f_matrix);

  rank = fmpz_mat_nullspace_mul(kernel, f_matrix_z);
  fmpz_mat_window_init(kernel_window, kernel, 0, 0, B.size(), rank);

  DEBUG_INFO(4,
    {
      printf("(-) kernel_window: ");
      fmpz_mat_print_dimensions(kernel_window);
      printf("\n");
    }
  )

  fmpz_mat_div_colwise_gcd(kernel_window);

  DEBUG_INFO(4,
    {
      printf("(-) kernel_window (cleared): ");
      fmpz_mat_print_dimensions(kernel_window);
      printf("\n");
    }
  )

  fmpz_mat_init(kernel_in_orig_basis, N_basis.size(), rank);
  fmpz_mat_mul(kernel_in_orig_basis, B_matrix_z, kernel_window);

  DEBUG_INFO(4,
    {
      printf("(-) kernel_in_orig_basis: ");
      fmpz_mat_print_dimensions(kernel_in_orig_basis);
      printf("\n");
    }
  )

  fmpz_mat_div_colwise_gcd(kernel_in_orig_basis);

  DEBUG_INFO(4,
    {
      printf("(-) kernel_in_orig_basis (cleared): ");
      fmpz_mat_print_dimensions(kernel_in_orig_basis);
      printf("\n");
    }
  )

  // TODO: we should just return the neg_space basis as a matrix instead
  // of this conversion to sparse rep.
  std::vector<ManinElement> neg_space;
  for (int col = 0; col < rank; col++) {
    std::vector<MBEWC> components;
    for (int row = 0; row < N_basis.size(); row++) {
      if (!(fmpz_is_zero(fmpz_mat_entry(kernel_in_orig_basis, row, col)))) {
        fmpq_t coeff;
        fmpq_init(coeff);
        fmpq_set_fmpz(coeff, fmpz_mat_entry(kernel_in_orig_basis, row, col));
        components.push_back(MBEWC(row, coeff));
        fmpq_clear(coeff);
      }
    }
    ManinElement element = ManinElement(N, components);
    element.mark_as_sorted_unchecked();
    neg_space.push_back(element);
  }

  DEBUG_INFO_PRINT(3, " neg_space dimension: %d\n", rank);

  fmpz_mat_window_clear(kernel_window);

  fmpq_mat_clear(neg_one);
  fmpq_mat_clear(f_matrix);

  fmpz_mat_clear(f_matrix_z);
  fmpz_mat_clear(kernel);
  fmpz_mat_clear(kernel_in_orig_basis);

  DEBUG_INFO_PRINT(2, "dim %zu -> +: %zu, -: %zu\n", B.size(), pos_space.size(), neg_space.size());

  return {.pos_space = pos_space, .neg_space = neg_space};
}

DecomposeResult DecomposeResult::empty() {
  return {
    .done = std::vector<std::vector<ManinElement>>(),
    .special = std::vector<std::vector<ManinElement>>(),
    .remaining = std::vector<std::vector<ManinElement>>()
  };
}

// TODO: should be faster to decompose using a matrix of the action on the newform subspace

DecomposeResult decompose(std::vector<ManinElement> B, FmpqMatrix& map_of_basis, bool dimension_only) {
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

  fmpq_mat_get_fmpz_mat_colwise(B_matrix_z, NULL, B_matrix);
  fmpq_mat_clear(B_matrix);

  DEBUG_INFO(4,
    {
      printf("B_matrix_z: ");
      fmpz_mat_print_dimensions(B_matrix_z);
      printf("\n");
    }
  )

  // Finds pivot rows of the matrix B.
  std::vector<int> pivots;
  fmpz* pivot_coeffs = _fmpz_vec_init(B.size());

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

  // Construct matrix of the linear map f acting on B
  fmpq_mat_t f_matrix;
  fmpq_mat_init(f_matrix, B.size(), B.size());

  fmpq_mat_t pivot_rows;
  fmpq_mat_init(pivot_rows, B.size(), N_basis.size());

  for (int row = 0; row < B.size(); row++) {
    for (int col = 0; col < N_basis.size(); col++) {
      fmpq_div_fmpz(
        fmpq_mat_entry(pivot_rows, row, col),
        fmpq_mat_entry(map_of_basis.mat, pivots[row], col),
        (pivot_coeffs + row)
      );
    }
  }

  DEBUG_INFO(4,
    {
      printf("pivot_rows: ");
      fmpq_mat_print_dimensions(pivot_rows);
      printf("\n");
    }
  )

  fmpq_mat_mul_fmpz_mat(f_matrix, pivot_rows, B_matrix_z);
  fmpq_mat_clear(pivot_rows);

  DEBUG_INFO(4,
    {
      printf("f_matrix: ");
      fmpq_mat_print_dimensions(f_matrix);
      printf("\n");
    }
  )

  // bool use_map_of_basis = true;

  // if (use_map_of_basis) {
  //   fmpq_mat_t map_of_basis;
  //   fmpq_mat_init(map_of_basis, B.size(), N_basis.size());
  //   for (int col = 0; col < N_basis.size(); col++) {
  //     ManinElement fb = f(N_basis[col]);

  //     printf("%zu out of %zu\n", fb.components.size(), N_basis.size());

  //     auto it = fb.components.begin();
  //     int pivot_index = 0;

  //     while (true) {
  //       if (it == fb.components.end()) break;
  //       if (pivot_index == pivots.size()) break;

  //       int row = it->basis_index;

  //       if (row > pivots[pivot_index]) {
  //         pivot_index++;
  //       } else if (row == pivots[pivot_index]) {
  //         fmpq_t coeff;
  //         fmpq_init(coeff);
  //         fmpq_set(coeff, it->coeff);
  //         fmpq_div_fmpz(coeff, coeff, (pivot_coeffs + pivot_index));
  //         fmpq_set(fmpq_mat_entry(map_of_basis, pivot_index, col), coeff);
  //         fmpq_clear(coeff);
  //         it++;
  //         pivot_index++;
  //       } else if (row < pivots[pivot_index]) {
  //         it++;
  //       }
  //     }
  //   }

  //   DEBUG_INFO(4,
  //     {
  //       printf("map_of_basis: ");
  //       fmpq_mat_print_dimensions(map_of_basis);
  //       printf("\n");
  //     }
  //   )

  //   fmpq_mat_mul_fmpz_mat(f_matrix, map_of_basis, B_matrix_z);
  //   fmpq_mat_clear(map_of_basis);

  //   DEBUG_INFO(4,
  //     {
  //       printf("f_matrix: ");
  //       fmpq_mat_print_dimensions(f_matrix);
  //       printf("\n");
  //     }
  //   )

  // } else {
  //   for (int col = 0; col < B.size(); col++) {
  //     ManinElement fb = B[col].map(f, N);

  //     // printf("rows: ");

  //     auto it = fb.components.begin();
  //     int pivot_index = 0;

  //     while (true) {
  //       if (it == fb.components.end()) break;
  //       if (pivot_index == pivots.size()) break;

  //       int row = it->basis_index;

  //       if (row > pivots[pivot_index]) {
  //         pivot_index++;
  //       } else if (row == pivots[pivot_index]) {
  //         fmpq_t coeff;
  //         fmpq_init(coeff);
  //         fmpq_set(coeff, it->coeff);
  //         fmpq_div_fmpz(coeff, coeff, (pivot_coeffs + pivot_index));
  //         fmpq_set(fmpq_mat_entry(f_matrix, pivot_index, col), coeff);
  //         fmpq_clear(coeff);
  //         it++;
  //         pivot_index++;
  //       } else if (row < pivots[pivot_index]) {
  //         it++;
  //       }
  //     }
  //   }
  // }

  _fmpz_vec_clear(pivot_coeffs, B.size());

  std::vector<std::vector<ManinElement>> done;
  std::vector<std::vector<ManinElement>> special;
  std::vector<std::vector<ManinElement>> remaining;

  fmpq_poly_t min_poly;
  fmpz_poly_t min_poly_z;
  fmpq_poly_init(min_poly);
  fmpz_poly_init(min_poly_z);


  if (fmpq_mat_is_zero(f_matrix)) {
    // XXX: FLINT seems to think that the minimal polynomial of the zero matrix is 1, and not T.
    // Manually set the minimal polynomial to x.
    fmpz_poly_set_coeff_si(min_poly_z, 1, 1);
  } else {
    fmpq_mat_minpoly(min_poly, f_matrix);
    fmpq_poly_get_numerator(min_poly_z, min_poly);
  }

  // Note: computing char_poly seems much slower than min_poly, so we're not using this code:
  // if (dimension_only) {
  //   if (fmpq_mat_is_zero(f_matrix)) {
  //     fmpz_poly_set_coeff_si(min_poly_z, B.size(), 1);
  //   } else {
  //     fmpq_mat_charpoly(min_poly, f_matrix);
  //     fmpq_poly_get_numerator(min_poly_z, min_poly);
  //   }

  int min_poly_degree = fmpz_poly_degree(min_poly_z);

  DEBUG_INFO_PRINT(3, " min_poly_z degree: %d\n", min_poly_degree);

  DEBUG_INFO(5,
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

  int num_factors = min_poly_factored->num;
  if (num_factors == 0) {
    // This should actually be impossible
    assert(false);
  } else if (num_factors == 1) {
    int deg = fmpz_poly_degree(min_poly_factored->p);
    if (deg == B.size()) {
      done.push_back(B);
      DEBUG_INFO_PRINT(3, " minimal polynomial irreducible and degree equal to space dimension %zu\n", B.size());
    } else {
      special.push_back(B);
      DEBUG_INFO_PRINT(3, " minimal polynomial irreducible but space dimension %zu is not equal to degree %d\n", B.size(), deg);
    }
  } else {
    fmpq_mat_t poly_on_f_matrix;
    fmpz_mat_t poly_on_f_matrix_z, poly_mat_kernel;
    fmpq_mat_init(poly_on_f_matrix, B.size(), B.size());
    fmpz_mat_init(poly_on_f_matrix_z, B.size(), B.size());
    fmpz_mat_init(poly_mat_kernel, B.size(), B.size());

    fmpz_mat_t poly_mat_kernel_window, poly_mat_kernel_in_orig_basis;

    int dimension_excess = 0;

    // For each factor g of the minpoly, we compute g(T) and take its kernel.
    for (int i = 0; i < num_factors; i++) {

      fmpz_poly_struct *factor = min_poly_factored->p + i;
      int exp = *(min_poly_factored->exp + i);
      int degree = fmpz_poly_degree(factor);

      if (dimension_only && (min_poly_degree + dimension_excess + degree > B.size())) {
        DEBUG_INFO_PRINT(3, " factor appears only once, degree: %d\n", degree);
        done.emplace_back(degree, ManinElement::zero(N));
        continue;
      }

      fmpz_poly_apply_fmpq_mat(poly_on_f_matrix, f_matrix, factor);

      DEBUG_INFO(4,
        {
          printf("poly_on_f_matrix: ");
          fmpq_mat_print_dimensions(poly_on_f_matrix);
          printf("\n");
        }
      )
      // NOTE: this needs to be rowwise!
      fmpq_mat_get_fmpz_mat_rowwise(poly_on_f_matrix_z, NULL, poly_on_f_matrix);

      DEBUG_INFO(4,
        {
          printf("poly_on_f_matrix_z: ");
          fmpz_mat_print_dimensions(poly_on_f_matrix_z);
          printf("\n");
        }
      )

      int rank = fmpz_mat_nullspace_mul(poly_mat_kernel, poly_on_f_matrix_z);

      fmpz_mat_window_init(poly_mat_kernel_window, poly_mat_kernel, 0, 0, B.size(), rank);

      DEBUG_INFO(4,
        {
          printf("poly_mat_kernel_window: ");
          fmpz_mat_print_dimensions(poly_mat_kernel_window);
          printf("\n");
        }
      )

      fmpz_mat_div_colwise_gcd(poly_mat_kernel_window);

      DEBUG_INFO(6,
        {
          fmpz_mat_print_pretty(poly_mat_kernel_window);
          printf("\n");
        }
      )

      DEBUG_INFO(4,
        {
          printf("poly_mat_kernel_window (cleared): ");
          fmpz_mat_print_dimensions(poly_mat_kernel_window);
          printf("\n");
        }
      )

      fmpz_mat_init(poly_mat_kernel_in_orig_basis, N_basis.size(), rank);
      fmpz_mat_mul(poly_mat_kernel_in_orig_basis, B_matrix_z, poly_mat_kernel_window);
      fmpz_mat_window_clear(poly_mat_kernel_window);

      DEBUG_INFO(4,
        {
          printf("poly_mat_kernel_in_orig_basis: ");
          fmpz_mat_print_dimensions(poly_mat_kernel_in_orig_basis);
          printf("\n");
        }
      )

      fmpz_mat_div_colwise_gcd(poly_mat_kernel_in_orig_basis);

      DEBUG_INFO(4,
        {
          printf("poly_mat_kernel_in_orig_basis (cleared): ");
          fmpz_mat_print_dimensions(poly_mat_kernel_in_orig_basis);
          printf("\n");
        }
      )

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
        dimension_excess += (rank - degree);
        remaining.push_back(output);
      }
    }

    fmpq_mat_clear(poly_on_f_matrix);
    fmpz_mat_clear(poly_on_f_matrix_z);
    fmpz_mat_clear(poly_mat_kernel);
  }

  fmpz_poly_clear(min_poly_z);
  fmpz_poly_factor_clear(min_poly_factored);

  fmpq_mat_clear(f_matrix);
  fmpz_mat_clear(B_matrix_z);

  return {.done = done, .special = special, .remaining = remaining};
}