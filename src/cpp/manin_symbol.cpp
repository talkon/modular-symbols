#include "manin_symbol.h"

#include "cache_decorator.h"

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fexpr_builtin.h>
#include <flint/fmpz_poly.h>
#include <flint/arith.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>

void ManinSymbol::print() {
  printf("(%lld, %lld)_%lld", c, d, N);
}

bool ManinSymbol::is_equivalent(const ManinSymbol& other) {
  if (N != other.N)
    return false;

  return ((c * other.d - d * other.c) % N) == 0;
}

// Base implementation of `manin_generators()`
// XXX: Surely there is a better approach to this pattern, but I think this works for now.
// TODO: Check that caching actually works across multiple compilation units.
std::vector<ManinGenerator> _impl_manin_generators(const int64_t n) {
  fmpz_t N;
  fmpz_init_set_ui(N, n);

  fmpz_poly_t divisors; // acts as a list of divisors of N
  fmpz_poly_init(divisors);
  arith_divisors(divisors, N);

  int64_t tau = fmpz_poly_length(divisors);

  int64_t index = 0;
  ManinGenerator first (index, {.N = n, .c = 0, .d = 1});
  index++;

  std::vector<ManinGenerator> out = {first};

  fmpz_t C, D, G, M, X;
  fmpz_init(C);
  fmpz_init(D);
  fmpz_init(G);
  fmpz_init(M);
  fmpz_init(X);

  for (int64_t i = 0; i < tau - 1; i++) {
    fmpz_poly_get_coeff_fmpz(D, divisors, i);
    fmpz_poly_get_coeff_fmpz(M, divisors, tau - 1 - i); // M = N / D
    int d = fmpz_get_si(D);
    int m = fmpz_get_si(M);
    // For each residue class x mod M, pick the smallest integer c such that gcd(D, c) = 1.
    // If gcd(D, M, x) > 1, we skip.
    for (int64_t x = 0; x < m; x++) {
      fmpz_set_ui(X, x);
      fmpz_gcd3(G, D, M, X);
      if (fmpz_cmp_si(G, 1) > 0) {
        continue;
      }
      for (int64_t c = x; c < n; c += m) {
        // printf("%d, %d\n", d, c);
        fmpz_set_ui(C, c);
        fmpz_gcd(G, D, C);
        if (fmpz_is_one(G)) {
          ManinGenerator new_gen (index, {.N = n, .c = d, .d = c});
          index++;
          // manin_generator_print(&new_gen);
          // printf("\n");
          out.push_back(new_gen);
          break;
        }
      }
    }

    for (int64_t x = 0; x < n; x++) {
      fmpz_set_ui(X, x);
      fmpz_gcd(G, D, X);
      if (fmpz_cmp_si(G, 1) > 0) {
        continue;
      }
    }
  }

  return out;
}

std::vector<ManinGenerator> manin_generators(const int64_t n) {
  static CacheDecorator<std::vector<ManinGenerator>, const int64_t> _cache_manin_generators(_impl_manin_generators);
  return _cache_manin_generators(n);
}

// Base implementation of `find_generator_index()`
ManinGenerator _impl_find_generator(const ManinSymbol ms) {
  std::vector<ManinGenerator> generators = manin_generators(ms.N);
  auto first = generators.begin();
  auto last = generators.end();

  auto mg = std::find_if(first, last,
    [&](ManinSymbol gen) { return gen.is_equivalent(ms); }
  );

  assert (mg != last); // Manin symbol should match one of the generators

  return *mg;
}

ManinGenerator find_generator(const ManinSymbol ms) {
  static CacheDecorator<ManinGenerator, const ManinSymbol> _cache_find_generator(_impl_find_generator);
  return _cache_find_generator(ms);
}

ManinGenerator ManinSymbol::as_generator() {
  return find_generator(*this);
}

int main() {
  std::vector<ManinGenerator> mgs;
  for (int i = 0; i < 10; i++) {
    mgs = manin_generators(2 * 3 * 5 * 7 * 11 * 13);
    printf("%zu\n", mgs.size());
  }
  mgs[2 * 3 * 5 * 7 * 11 * 13 + 50000].print();
  printf("\n");

  for (int i = 11; i < 20; i++) {
    ManinSymbol ms1 = {.N = i, .c = 4, .d = 7};
    ManinGenerator mg1 = find_generator(ms1);
    ms1.print();
    mg1.print();
    printf(" %lld %p\n", mg1.index, &mg1);
  }

  return 0;
}