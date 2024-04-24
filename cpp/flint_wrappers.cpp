#include "flint_wrappers.h"

void FmpqMatrix::set_copy(const fmpq_mat_t src) {
  fmpq_mat_init_set(mat, src);
}

void FmpqMatrix::set_move(fmpq_mat_t src) {
  fmpq_mat_init(mat, 0, 0);
  fmpq_mat_swap(mat, src);
  fmpq_mat_clear(src);
}

FmpqMatrix::FmpqMatrix(const FmpqMatrix& src) {
  fmpq_mat_init_set(mat, src.mat);
}

FmpqMatrix::FmpqMatrix(FmpqMatrix&& src) {
  fmpq_mat_init(mat, 0, 0);
  fmpq_mat_swap(mat, src.mat);
}

FmpqMatrix& FmpqMatrix::operator=(const FmpqMatrix& src) {
  fmpq_mat_init_set(mat, src.mat);
  return *this;
}

FmpqMatrix& FmpqMatrix::operator=(FmpqMatrix&& src) {
  fmpq_mat_init(mat, 0, 0);
  fmpq_mat_swap(mat, src.mat);
  return *this;
}

FmpqMatrix::~FmpqMatrix() {
  fmpq_mat_clear(mat);
}

void FmpzPoly::set_copy(const fmpz_poly_t src) {
  fmpz_poly_init(poly);
  fmpz_poly_set(poly, src);
}

void FmpzPoly::set_move(fmpz_poly_t src) {
  fmpz_poly_init(poly);
  fmpz_poly_swap(poly, src);
  fmpz_poly_clear(src);
}

FmpzPoly::FmpzPoly(const FmpzPoly& src) {
  fmpz_poly_init(poly);
  fmpz_poly_set(poly, src.poly);
}

FmpzPoly::FmpzPoly(FmpzPoly&& src) {
  fmpz_poly_init(poly);
  fmpz_poly_swap(poly, src.poly);
}

FmpzPoly& FmpzPoly::operator=(const FmpzPoly& src) {
  fmpz_poly_init(poly);
  fmpz_poly_set(poly, src.poly);
  return *this;
}

FmpzPoly& FmpzPoly::operator=(FmpzPoly&& src) {
  fmpz_poly_init(poly);
  fmpz_poly_swap(poly, src.poly);
  return *this;
}

FmpzPoly::~FmpzPoly() {
  fmpz_poly_clear(poly);
}
