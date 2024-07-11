#ifndef FLINT_WRAPPERS_H
#define FLINT_WRAPPERS_H

#include <flint/fmpq_mat.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>

// C++ wrappers for FLINT datatypes, so that we can easily pass around FLINT matrices and put them in C++ vectors, etc.

struct FmpqMatrix {
  fmpq_mat_t mat;
  bool initialized = false;

  FmpqMatrix() {};
  void set_copy(const fmpq_mat_t);
  void set_move(fmpq_mat_t);

  // Copy and move constructors/operators
  FmpqMatrix(const FmpqMatrix&);
  FmpqMatrix(FmpqMatrix&&);
  FmpqMatrix& operator=(const FmpqMatrix&);
  FmpqMatrix& operator=(FmpqMatrix&&);

  void clear();
  ~FmpqMatrix();
};

struct FmpzMatrix {
  fmpz_mat_t mat;
  bool initialized = false;

  FmpzMatrix() {};
  void set_copy(const fmpz_mat_t);
  void set_move(fmpz_mat_t);

  // Copy and move constructors/operators
  FmpzMatrix(const FmpzMatrix&);
  FmpzMatrix(FmpzMatrix&&);
  FmpzMatrix& operator=(const FmpzMatrix&);
  FmpzMatrix& operator=(FmpzMatrix&&);

  void clear();
  ~FmpzMatrix();
};


struct FmpzPoly {
  fmpz_poly_t poly;
  bool initialized = false;

  FmpzPoly() {};
  void set_copy(const fmpz_poly_t);
  void set_move(fmpz_poly_t);

  // Copy and move constructors/operators
  FmpzPoly(const FmpzPoly&);
  FmpzPoly(FmpzPoly&&);
  FmpzPoly& operator=(const FmpzPoly&);
  FmpzPoly& operator=(FmpzPoly&&);

  void clear();
  ~FmpzPoly();
};

#endif // FLINT_WRAPPERS_H