#include "modular_symbol.h"
#include "manin_element.h"
#include "manin_basis.h"
#include "newspace.h"
#include "linalg.h"
#include "utils.h"
#include "debug_utils.h"
#include "cache_decorator.h"

#include <flint/ulong_extras.h>
#include <cmath>

std::vector<IntMatrix2x2> heilbronn_matrices(int64_t p) {
  std::vector<IntMatrix2x2> result;
  result.push_back({.x = 1, .y = 0, .z = 0, .w = p});

  for (int r = 0; r < p; r++) {
    int x1 = p;
    int x2 = -r;
    int y1 = 0;
    int y2 = 1;
    int a = -p;
    int b = r;
    result.push_back({.x = x1, .y = x2, .z = y1, .w = y2});
    while (b != 0) {
      // XXX: This rounding behavior is slightly inefficient (gives slightly more matrices than necessary), but is correct afaict.
      int q = std::lround((double) a / (double) b);
      int c = a - b * q;
      a = -b;
      b = c;
      int x3 = q * x2 - x1;
      x1 = x2;
      x2 = x3;
      int y3 = q * y2 - y1;
      y1 = y2;
      y2 = y3;
      result.push_back({.x = x1, .y = x2, .z = y1, .w = y2});
    }
  }

  return result;
}

// TODO: consider caching this, or really, just caching its matrix.

ManinElement _impl_hecke_action(ManinBasisElement mbe, int64_t p) {
  int64_t level = mbe.N;
  ManinElement result = ManinElement::zero(level);
  for (auto mat : heilbronn_matrices(p)) {
    ManinGenerator mg = mbe.right_action_by(mat).as_generator();
    // printf("mg: "); mg.print();
    // printf("\n");
    // auto mbe = level_and_index_to_basis(level, mg.index);
    // printf("mbe: "); mbe.print_with_generators();
    // printf("\n");
    result += level_and_index_to_basis(level, mg.index);
  }
  // printf("result: "); result.print_with_generators();
  // printf("\n\n");
  return result;
}

ManinElement hecke_action(ManinBasisElement mbe, int64_t p) {
  static CacheDecorator<ManinElement, ManinBasisElement, int64_t> _cache_hecke_action(_impl_hecke_action);
  return _cache_hecke_action(mbe, p);
}

ManinElement _impl_atkin_lehner_action(ManinBasisElement mbe, int64_t q) {
  int64_t level = mbe.N;
  utils::XgcdResult xgcd = utils::xgcd(q, level / q);
  IntMatrix2x2 matrix = {.x = q, .y = -xgcd.b, .z = level, .w = xgcd.a * q};
  return mbe.as_modular_symbol().left_action_by(matrix).to_manin_element(level);
}

ManinElement atkin_lehner_action(ManinBasisElement mbe, int64_t q) {
  static CacheDecorator<ManinElement, ManinBasisElement, int64_t> _cache_atkin_lehner_action(_impl_atkin_lehner_action);
  return _cache_atkin_lehner_action(mbe, q);
}

std::vector<std::vector<ManinElement>> newform_subspaces(int64_t level, bool use_atkin_lehner) {
  std::vector<ManinElement> basis = newspace_basis(level);

  DEBUG_INFO_PRINT(1, "Starting computation of newform subspaces for level %lld\n", level);

  std::vector<std::vector<ManinElement>> done;
  std::vector<std::vector<ManinElement>> remaining = { basis };

  if (use_atkin_lehner) {
    n_factor_t factors;
    n_factor_init(&factors);
    n_factor(&factors, level, 1);
    for (int i = 0; i < factors.num && remaining.size() > 0; i++){
      int64_t q = n_pow(factors.p[i], factors.exp[i]);
      DEBUG_INFO_PRINT(1, "Decomposing spaces using Atkin-Lehner involution for prime %lld\n", (int64_t) factors.p[i]);
      auto f = [q](ManinBasisElement mbe) { return atkin_lehner_action(mbe, q); };
      std::vector<std::vector<ManinElement>> new_remaining;
      // XXX: This causes the action of `f` to be recomputed many times.
      for (auto subspace_basis : remaining) {
        DEBUG_INFO_PRINT(2, "Space size: %zu\n", subspace_basis.size());
        // TODO: since all eigenvalues are +1 or -1, we don't need to use the full power of decompose()
        DecomposeResult dr = decompose(subspace_basis, f);
        // minpoly factor degree should always be 1, so if anything goes in done, it's actually done
        done.insert(done.end(), dr.done.begin(), dr.done.end());
        new_remaining.insert(new_remaining.end(), dr.remaining.begin(), dr.remaining.end());
      }
      remaining = new_remaining;
    }
  }

  // n_factor is just a struct, does not need to be cleared

  n_primes_t prime_iter;
  n_primes_init(prime_iter);
  while (remaining.size() > 0) {
    int64_t p = n_primes_next(prime_iter);
    if (level % p == 0) continue;

    DEBUG_INFO_PRINT(1, "decomposing spaces using Hecke operator for prime %lld\n", p);
    auto f = [p](ManinBasisElement mbe) { return hecke_action(mbe, p); };
    std::vector<std::vector<ManinElement>> new_remaining;
    // XXX: This causes the action of `f` to be recomputed many times.
    for (auto subspace_basis : remaining) {
      DEBUG_INFO_PRINT(2, "space size: %zu\n", subspace_basis.size());
      DecomposeResult dr = decompose(subspace_basis, f);
      done.insert(done.end(), dr.done.begin(), dr.done.end());
      new_remaining.insert(new_remaining.end(), dr.remaining.begin(), dr.remaining.end());
    }
    remaining = new_remaining;
  }

  n_primes_clear(prime_iter);

  return done;
}

std::vector<int> newform_subspace_dimensions(int64_t level, bool use_atkin_lehner) {
  auto nss = newform_subspaces(level, use_atkin_lehner);
  std::vector<int> sizes;
  for (auto ns : nss) {
    sizes.push_back(ns.size());
  }
  std::sort(sizes.begin(), sizes.end());
  return sizes;
}
