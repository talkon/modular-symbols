#include "linalg.h"
#include "manin_basis.h"
#include "manin_element.h"

#include <flint/fmpq_mat.h>
#include <flint/fmpz_mat.h>

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

  // XXX: this feels a bit wasteful
  fmpq_mat_get_fmpz_mat_colwise(B_matrix_z, NULL, B_matrix);
  fmpq_mat_clear(B_matrix);

  // printf("basis matrix:\n");
  // fmpz_mat_print_pretty(B_matrix_z);
  // printf("\n");

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
  fmpq_mat_get_fmpz_mat_rowwise(map_matrix_z, NULL, map_matrix);
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

  // FIXME: multiply result matrix with basis B first.

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