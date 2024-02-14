#ifndef UTILS_H
#define UTILS_H

#include <flint/fmpz.h>

#include <cstdint>

struct XgcdResult {
  int64_t d;
  int64_t a;
  int64_t b;
};

// Computes a, b, d such that d = gcd(x, y) and d = ax + by.
inline XgcdResult xgcd(int64_t x, int64_t y) {
  fmpz_t A, B, X, Y, D;

  fmpz_init(A);
  fmpz_init(B);
  fmpz_init(D);
  fmpz_init_set_si(X, x);
  fmpz_init_set_si(Y, y);

  fmpz_xgcd_canonical_bezout(D, A, B, X, Y);

  XgcdResult result = {
    .d = fmpz_get_si(D),
    .a = fmpz_get_si(A),
    .b = fmpz_get_si(B)
  };

  fmpz_clear(A);
  fmpz_clear(B);
  fmpz_clear(D);
  fmpz_clear(X);
  fmpz_clear(Y);

  return result;
}

#endif // UTILS_H