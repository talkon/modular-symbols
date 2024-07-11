#include "subspace_basis.h"

#include "manin_basis.h"

#include <cassert>

DenseBasis sparse_to_dense(SparseBasis& sparse, int64_t level, bool clear) {

  int nrows = manin_basis(level).size();
  int ncols = sparse.size();

  fmpq_mat_t dense_col_q;
  fmpz_mat_t dense_mat_z, dense_mat_z_window;

  fmpq_mat_init(dense_col_q, nrows, 1);
  fmpz_mat_init(dense_mat_z, nrows, ncols);
  fmpz_mat_zero(dense_mat_z);

  for (int col = 0; col < ncols; col++) {
    fmpq_mat_zero(dense_col_q);
    for (auto& component : sparse[col].components) {
      int row = component.basis_index;
      assert(row < nrows);
      fmpq_set(fmpq_mat_entry(dense_col_q, row, 0), component.coeff);
    }

    fmpz_mat_window_init(dense_mat_z_window, dense_mat_z, 0, col, nrows, col + 1);
    fmpq_mat_get_fmpz_mat_colwise(dense_mat_z_window, NULL, dense_col_q);
    fmpz_mat_window_clear(dense_mat_z_window);
  }

  fmpq_mat_clear(dense_col_q);

  if (clear) sparse.clear();

  DenseBasis dense;
  dense.set_move(dense_mat_z);
  return dense;
}

SparseBasis dense_to_sparse(DenseBasis& dense, int64_t level, bool clear) {

  int nrows = dense.mat->r;
  int ncols = dense.mat->c;

  SparseBasis sparse;

  for (int col = 0; col < ncols; col++) {
    std::vector<MBEWC> components;
    for (int row = 0; row < nrows; row++) {
      if (!(fmpz_is_zero(fmpz_mat_entry(dense.mat, row, col)))) {
        fmpq_t coeff;
        fmpq_init(coeff);
        fmpq_set_fmpz(coeff, fmpz_mat_entry(dense.mat, row, col));
        components.emplace_back(row, coeff);
        fmpq_clear(coeff);
      }
    }
    // has_duplicate_keys(components);
    ManinElement element = ManinElement(level, components);
    element.mark_as_sorted_unchecked();
    sparse.push_back(element);
  }

  if (clear) dense.clear();

  return sparse;
}