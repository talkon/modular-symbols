#ifndef FLINT_WRAPPERS_H
#define FLINT_WRAPPERS_H

#include <flint/fmpq_mat.h>
#include <flint/fmpz_poly.h>

// C++ wrappers for flint datatypes

struct FmpqMatrix {
  fmpq_mat_t mat;

  FmpqMatrix() {};
  void set_copy(const fmpq_mat_t);
  void set_move(fmpq_mat_t);

  // Copy and move constructors/operators
  FmpqMatrix(const FmpqMatrix&);
  FmpqMatrix(FmpqMatrix&&);
  FmpqMatrix& operator=(const FmpqMatrix&);
  FmpqMatrix& operator=(FmpqMatrix&&);

  ~FmpqMatrix();
};

struct FmpzPoly {
  fmpz_poly_t poly;

  FmpzPoly() {};
  void set_copy(const fmpz_poly_t);
  void set_move(fmpz_poly_t);

  // Copy and move constructors/operators
  FmpzPoly(const FmpzPoly&);
  FmpzPoly(FmpzPoly&&);
  FmpzPoly& operator=(const FmpzPoly&);
  FmpzPoly& operator=(FmpzPoly&&);

  ~FmpzPoly();
};

#endif // FLINT_WRAPPERS_H