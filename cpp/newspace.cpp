#include "newspace.h"
#include "boundary_map.h"
#include "modular_symbol.h"
#include "manin_symbol.h"
#include "manin_element.h"
#include "linalg.h"
#include "debug_utils.h"
#include "cache_decorator.h"

#include <flint/fmpz_poly.h>
#include <flint/ulong_extras.h>

#include <cassert>

// #define CACHE_OLDSPACE_MAP

#ifdef CACHE_OLDSPACE_MAP
ManinElement _impl_oldspace_map(ManinBasisElement mbe, int64_t d, int64_t M) {
#else
ManinElement oldspace_map(ManinBasisElement mbe, int64_t d, int64_t M) {
#endif
  int64_t N = mbe.N;
  assert(N % (d * M) == 0);
  IntMatrix2x2 matrix = {.x = d, .y = 0, .z = 0, .w = 1};
  return mbe.as_modular_symbol().left_action_by(matrix).to_manin_element(M);
}

#ifdef CACHE_OLDSPACE_MAP
ManinElement oldspace_map(ManinBasisElement mbe, int64_t d, int64_t M) {
  static CacheDecorator<ManinElement, ManinBasisElement, int64_t, int64_t> _cache_oldspace_map(_impl_oldspace_map);
  return _cache_oldspace_map(mbe, d, M);
}
#endif

std::vector<ManinElement> newspace_basis(int64_t level) {
  std::vector<ManinElement> current_basis = cuspidal_manin_basis(level);

  DEBUG_INFO_PRINT(1, "Started computation of newspace basis for level %lld\n", level)
  DEBUG_INFO_PRINT(2, "Starting current_basis size: %zu\n", current_basis.size())

  n_factor_t factors;
  n_factor_init(&factors);
  n_factor(&factors, level, 1);

  // Seems like it suffices to compute the oldspace maps only for M of the form N / p, where p is a prime factor of N.

  // Threshold to decide for each m whether to first compute d = 1 then d = p in two separate passes or compute both
  // at the same time. It seems like this should just be 0 (i.e. for everything, always compute d = 1 maps before d = p maps).
  int threshold = 0;

  // In general, p is sorted by size, so that m is sorted large to small, so that `current_basis` gets smaller faster.

  // We first compute all the d = 1 oldspace maps for these M because the d = 1 maps lead to less
  // blowup in the size of matrix entries
  for (int i = 0; i < factors.num; i++) {
    int64_t p = factors.p[i];
    int64_t m = level / p;

    if (m < threshold || m < 11 || m == 12 || m == 13 || m == 16 || m == 18 || m == 25 ) {
      continue;
    }

    DEBUG_INFO_PRINT(2, "Computing d = 1 oldspace map for M: %lld\n", m);

    auto f1 = [m](ManinBasisElement mbe) { return oldspace_map(mbe, 1, m); };
    current_basis = map_kernel(current_basis, f1, m);
    DEBUG_INFO_PRINT(3, "M: %lld, d: 1, current_basis size: %zu\n", m, current_basis.size());
    DEBUG_INFO_PRINT(4, "current_basis alloc_size: %lu bytes\n", ManinElement::vector_alloc_size(current_basis));
  }

  // Compute the d = p oldspace maps for large M
  for (int i = 0; i < factors.num; i++) {
    int64_t p = factors.p[i];
    int64_t m = level / p;

    if (m < threshold || m < 11 || m == 12 || m == 13 || m == 16 || m == 18 || m == 25 ) {
      continue;
    }

    DEBUG_INFO_PRINT(2, "Computing d = p oldspace map for M: %lld\n", m);

    auto fp = [p, m](ManinBasisElement mbe) { return oldspace_map(mbe, p, m); };
    current_basis = map_kernel(current_basis, fp, m);
    DEBUG_INFO_PRINT(3, "M: %lld, d: %lld, current_basis size: %zu\n", m, p, current_basis.size());
    DEBUG_INFO_PRINT(4, "current_basis alloc_size: %lu bytes\n", ManinElement::vector_alloc_size(current_basis));
  }

  // Compute maps for "small" M, i.e. M < threshold (This does not seem to help)
  // for (int i = 0; i < factors.num; i++) {
  //   int64_t p = factors.p[i];
  //   int64_t m = level / p;

  //   // Skip values of m such that m and all its divisors have trivial newspaces
  //   if (m > threshold || m < 11 || m == 12 || m == 13 || m == 16 || m == 18 || m == 25 ) {
  //     continue;
  //   }

  //   DEBUG_INFO_PRINT(2, "Computing oldspace maps for M: %lld\n", m);

  //   auto f1 = [m](ManinBasisElement mbe) { return oldspace_map(mbe, 1, m); };
  //   current_basis = map_kernel(current_basis, f1, m);
  //   DEBUG_INFO_PRINT(3, "M: %lld, d: 1, current_basis size: %zu\n", m, current_basis.size());

  //   auto fp = [p, m](ManinBasisElement mbe) { return oldspace_map(mbe, p, m); };
  //   current_basis = map_kernel(current_basis, fp, m);
  //   DEBUG_INFO_PRINT(3, "M: %lld, d: %lld, current_basis size: %zu\n", m, p, current_basis.size());
  // }

  return current_basis;
}