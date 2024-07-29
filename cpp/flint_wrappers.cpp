#include "flint_wrappers.h"

void FmpqMatrix::set_copy(const fmpq_mat_t src) {
  fmpq_mat_init_set(mat, src);
  initialized = true;
}

void FmpqMatrix::set_move(fmpq_mat_t src) {
  fmpq_mat_init(mat, 0, 0);
  fmpq_mat_swap(mat, src);
  fmpq_mat_clear(src);
  initialized = true;
}

FmpqMatrix::FmpqMatrix(const FmpqMatrix& src) {
  fmpq_mat_init_set(mat, src.mat);
  initialized = true;
}

FmpqMatrix::FmpqMatrix(FmpqMatrix&& src) {
  fmpq_mat_init(mat, 0, 0);
  fmpq_mat_swap(mat, src.mat);
  initialized = true;
}

FmpqMatrix& FmpqMatrix::operator=(const FmpqMatrix& src) {
  fmpq_mat_init_set(mat, src.mat);
  initialized = true;
  return *this;
}

FmpqMatrix& FmpqMatrix::operator=(FmpqMatrix&& src) {
  fmpq_mat_init(mat, 0, 0);
  fmpq_mat_swap(mat, src.mat);
  initialized = true;
  return *this;
}

void FmpqMatrix::clear() {
  if (initialized) {
    fmpq_mat_clear(mat);
    initialized = false;
  }
}

FmpqMatrix::~FmpqMatrix() {
  if (initialized)
    fmpq_mat_clear(mat);
}

void FmpzMatrix::set_copy(const fmpz_mat_t src) {
  fmpz_mat_init_set(mat, src);
  initialized = true;
}

void FmpzMatrix::set_move(fmpz_mat_t src) {
  fmpz_mat_init(mat, 0, 0);
  fmpz_mat_swap(mat, src);
  fmpz_mat_clear(src);
  initialized = true;
}

FmpzMatrix::FmpzMatrix(const FmpzMatrix& src) {
  fmpz_mat_init_set(mat, src.mat);
  initialized = true;
}

FmpzMatrix::FmpzMatrix(FmpzMatrix&& src) {
  fmpz_mat_init(mat, 0, 0);
  fmpz_mat_swap(mat, src.mat);
  initialized = true;
}

FmpzMatrix& FmpzMatrix::operator=(const FmpzMatrix& src) {
  fmpz_mat_init_set(mat, src.mat);
  initialized = true;
  return *this;
}

FmpzMatrix& FmpzMatrix::operator=(FmpzMatrix&& src) {
  fmpz_mat_init(mat, 0, 0);
  fmpz_mat_swap(mat, src.mat);
  initialized = true;
  return *this;
}

void FmpzMatrix::clear() {
  if (initialized) {
    fmpz_mat_clear(mat);
    initialized = false;
  }
}

FmpzMatrix::~FmpzMatrix() {
  if (initialized)
    fmpz_mat_clear(mat);
}

void FmpzPoly::set_copy(const fmpz_poly_t src) {
  fmpz_poly_init(poly);
  fmpz_poly_set(poly, src);
  initialized = true;
}

void FmpzPoly::set_move(fmpz_poly_t src) {
  fmpz_poly_init(poly);
  fmpz_poly_swap(poly, src);
  fmpz_poly_clear(src);
  initialized = true;
}

FmpzPoly::FmpzPoly(const FmpzPoly& src) {
  fmpz_poly_init(poly);
  fmpz_poly_set(poly, src.poly);
  initialized = true;
}

FmpzPoly::FmpzPoly(FmpzPoly&& src) {
  fmpz_poly_init(poly);
  fmpz_poly_swap(poly, src.poly);
  initialized = true;
}

FmpzPoly& FmpzPoly::operator=(const FmpzPoly& src) {
  fmpz_poly_init(poly);
  fmpz_poly_set(poly, src.poly);
  initialized = true;
  return *this;
}

FmpzPoly& FmpzPoly::operator=(FmpzPoly&& src) {
  fmpz_poly_init(poly);
  fmpz_poly_swap(poly, src.poly);
  initialized = true;
  return *this;
}

void FmpzPoly::clear() {
  if (initialized) {
    fmpz_poly_clear(poly);
    initialized = false;
  }
}

FmpzPoly::~FmpzPoly() {
  if (initialized)
    fmpz_poly_clear(poly);
}