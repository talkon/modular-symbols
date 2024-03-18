#include "newspace.h"
#include "boundary_map.h"
#include "modular_symbol.h"
#include "manin_symbol.h"
#include "manin_element.h"
#include "linalg.h"
#include "debug_utils.h"
#include "cache_decorator.h"

#include <flint/fmpz_poly.h>
#include <flint/arith.h>

#include <cassert>

ManinElement oldspace_map(ManinBasisElement mbe, int64_t d, int64_t M) {
  int64_t N = mbe.N;
  assert(N % (d * M) == 0);
  IntMatrix2x2 matrix = {.x = d, .y = 0, .z = 0, .w = 1};
  return mbe.as_modular_symbol().left_action_by(matrix).to_manin_element(M);
}

// ManinElement oldspace_map(ManinBasisElement mbe, int64_t d, int64_t M) {
//   static CacheDecorator<ManinElement, ManinBasisElement, int64_t, int64_t> _cache_oldspace_map(_impl_oldspace_map);
//   return _cache_oldspace_map(mbe, d, M);
// }

std::vector<ManinElement> newspace_basis(int64_t level) {
  std::vector<ManinElement> current_basis = cuspidal_manin_basis(level);

  DEBUG_INFO_PRINT(1, "Started computation of newspace basis for level %lld\n", level)
  DEBUG_INFO_PRINT(2, "Starting current_basis size: %zu\n", current_basis.size())

  fmpz_t N, M, D, N_over_M;
  fmpz_init_set_si(N, level);
  fmpz_init(M);
  fmpz_init(D);
  fmpz_init(N_over_M);

  fmpz_poly_t divisors; // acts as a list of divisors of N
  fmpz_poly_t divisors_N_over_M; // acts as a list of divisors of N/M
  fmpz_poly_init(divisors);
  fmpz_poly_init(divisors_N_over_M);
  arith_divisors(divisors, N);

  int64_t tau = fmpz_poly_length(divisors);
  for (int i = 0; i < tau - 1; i++) {
    fmpz_poly_get_coeff_fmpz(M, divisors, i);
    int64_t m = fmpz_get_si(M);
    if (m < 11) {
      continue;
    }

    fmpz_poly_get_coeff_fmpz(N_over_M, divisors, tau - i - 1);
    arith_divisors(divisors_N_over_M, N_over_M);
    int64_t tau_N_over_M = fmpz_poly_length(divisors_N_over_M);
    DEBUG_INFO_PRINT(2, "Computing oldspace maps for M: %lld\n", m);
    for (int j = 0; j < tau_N_over_M; j++) {
      fmpz_poly_get_coeff_fmpz(D, divisors_N_over_M, j);
      int64_t d = fmpz_get_si(D);
      auto f = [d, m](ManinBasisElement mbe) { return oldspace_map(mbe, d, m); };
      current_basis = map_kernel(current_basis, f, m);
      DEBUG_INFO_PRINT(3, "M: %lld, d: %lld, current_basis size: %zu\n", m, d, current_basis.size());
    }
  }

  fmpz_clear(N);
  fmpz_clear(M);
  fmpz_clear(D);
  fmpz_clear(N_over_M);
  fmpz_poly_clear(divisors);
  fmpz_poly_clear(divisors_N_over_M);

  return current_basis;
}