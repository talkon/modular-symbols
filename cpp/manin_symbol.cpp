#include "manin_symbol.h"
#include "manin_basis.h"
#include "manin_element.h"
#include "modular_symbol.h"
#include "utils.h"

#include "cache_decorator.h"

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>
#include <flint/fexpr_builtin.h>
#include <flint/arith.h>

#include <iostream>
#include <iterator>
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

ManinSymbol ManinSymbol::apply_eta() {
  return {.N = this->N, .c = -(this->c), .d = this->d};
}

ManinSymbol ManinSymbol::apply_S() {
  return {.N = this->N, .c = this->d, .d = this->c};
}

ManinSymbol ManinSymbol::apply_T() {
  return {.N = this->N, .c = this->c + this->d, .d = -(this->c)};
}

ManinSymbol ManinSymbol::apply_T_2() {
  return {.N = this->N, .c = this->d, .d = -(this->c + this->d)};
}

ModularSymbol ManinSymbol::as_modular_symbol() {
  utils::XgcdResult xgcd = utils::xgcd(c, d);
  return {.a = xgcd.b, .b = -xgcd.a, .c = c, .d = d};
}

ManinElement ManinGenerator::as_element_unchecked() {
  fmpq_t one;
  fmpq_init(one);
  fmpq_one(one);
  std::vector<MGWC> components = {{.coeff = *one, .index = index}};

  ManinElement result = {.N = N, .components = components};
  result.mark_as_sorted_unchecked();
  return result;
}


// Base implementation of `manin_generators()`
// XXX: Surely there is a better approach to this caching pattern, but I think this works for now.
// XXX: It seems like this caching approach works across multiple compilation units,
// but I'm not sure why.
std::vector<ManinGenerator> _impl_manin_generators(const int64_t level) {
  fmpz_t N;
  fmpz_init_set_ui(N, level);

  fmpz_poly_t divisors; // acts as a list of divisors of N
  fmpz_poly_init(divisors);
  arith_divisors(divisors, N);

  int64_t tau = fmpz_poly_length(divisors);

  int64_t index = 0;
  ManinGenerator first (index, {.N = level, .c = 0, .d = 1});
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
      for (int64_t c = x; c < level; c += m) {
        // printf("%d, %d\n", d, c);
        fmpz_set_ui(C, c);
        fmpz_gcd(G, D, C);
        if (fmpz_is_one(G)) {
          ManinGenerator new_gen (index, {.N = level, .c = d, .d = c});
          index++;
          // manin_generator_print(&new_gen);
          // printf("\n");
          out.push_back(new_gen);
          break;
        }
      }
    }

    for (int64_t x = 0; x < level; x++) {
      fmpz_set_ui(X, x);
      fmpz_gcd(G, D, X);
      if (fmpz_cmp_si(G, 1) > 0) {
        continue;
      }
    }
  }

  fmpz_clear(C);
  fmpz_clear(D);
  fmpz_clear(G);
  fmpz_clear(M);
  fmpz_clear(X);
  fmpz_poly_clear(divisors);

  return out;
}

std::vector<ManinGenerator> manin_generators(const int64_t level) {
  // [ ]: figure out a way to share this variable across translation units?
  static CacheDecorator<std::vector<ManinGenerator>, const int64_t> _cache_manin_generators(_impl_manin_generators);
  return _cache_manin_generators(level);
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