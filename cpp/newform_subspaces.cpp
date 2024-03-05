#include "modular_symbol.h"
#include "manin_element.h"
#include "manin_basis.h"
#include "newspace.h"
#include "linalg.h"
#include "debug_timer.h"

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
ManinElement hecke_action(ManinBasisElement mbe, int64_t p) {
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

std::vector<std::vector<ManinElement>> newform_subspaces(int64_t level) {
  std::vector<ManinElement> basis = newspace_basis(level);

  info_with_time();
  printf(" starting computation of newform subspaces for level %lld\n", level);
  // printf("basis:\n");
  // for(auto elt : basis) {
  //   elt.print();
  //   printf("\n\n");
  // }

  std::vector<std::vector<ManinElement>> done;
  std::vector<std::vector<ManinElement>> remaining = { basis };

  // TODO: add Atkin-Lehner

  n_primes_t prime_iter;
  n_primes_init(prime_iter);
  while (remaining.size() > 0) {
    int64_t p = n_primes_next(prime_iter);
    if (level % p == 0) continue;

    info_with_time();
    printf(" decomposing spaces using prime %lld\n", p);
    auto f = [p](ManinBasisElement mbe) { return hecke_action(mbe, p); };
    std::vector<std::vector<ManinElement>> new_remaining;
    // XXX: This causes the action of `f` to be recomputed many times.
    for (auto subspace_basis : remaining) {
      info_with_time();
      printf(" space size: %zu\n", subspace_basis.size());
      DecomposeResult dr = decompose(subspace_basis, f);
      done.insert(done.end(), dr.done.begin(), dr.done.end());
      new_remaining.insert(new_remaining.end(), dr.remaining.begin(), dr.remaining.end());
    }
    remaining = new_remaining;
  }

  n_primes_clear(prime_iter);

  return done;
}

std::vector<int> newform_subspace_dimensions(int64_t level) {
  auto nss = newform_subspaces(level);
  std::vector<int> sizes;
  for (auto ns : nss) {
    sizes.push_back(ns.size());
  }
  std::sort(sizes.begin(), sizes.end());
  return sizes;
}
