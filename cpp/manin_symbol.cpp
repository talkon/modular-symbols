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

void ManinSymbol::print() const {
  printf("(%lld, %lld)_%lld", c, d, N);
}

bool ManinSymbol::operator<(const ManinSymbol& other) const {
  return (std::tuple<int64_t,int64_t,int64_t>) {N, c, d} < (std::tuple<int64_t,int64_t,int64_t>) {other.N, other.c, other.d};
}

bool ManinSymbol::is_equivalent(const ManinSymbol& other) const {
  if (N != other.N)
    return false;

  return ((c * other.d - d * other.c) % N) == 0;
}

ManinSymbol ManinSymbol::repr() {
  return {.N = this->N, .c = ((this->c) % N + N) % N, .d = ((this->d) % N + N) % N};
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

ManinSymbol ManinSymbol::right_action_by(IntMatrix2x2 mat) {
  return {.N = this->N, .c = this->c * mat.x + this->d * mat.z, .d = this->c * mat.y + this->d * mat.w};
}

ModularSymbol ManinSymbol::as_modular_symbol() {
  utils::XgcdResult xgcd = utils::xgcd(c, d);
  return {.a = xgcd.b, .b = -xgcd.a, .c = c, .d = d};
}

struct GeneratorComputationResult {
  std::vector<ManinGenerator> generators;
  std::map<ManinSymbol, int> generator_to_index;
  std::map<std::pair<int, int>, int> c_and_mod_M_to_d;

  GeneratorComputationResult(
    std::vector<ManinGenerator> generators,
    std::map<ManinSymbol, int> generator_to_index,
    std::map<std::pair<int, int>, int> c_and_mod_M_to_d
  ) :
    generators(generators),
    generator_to_index(generator_to_index),
    c_and_mod_M_to_d(c_and_mod_M_to_d)
  {}
};

// Base implementation of `manin_generators()`
// XXX: Surely there is a better approach to this caching pattern, but I think this works for now.
// XXX: It seems like this caching approach works across multiple compilation units,
// but I'm not sure why.
GeneratorComputationResult _impl_compute_manin_generators(const int64_t level) {
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
  std::map<std::pair<int, int>, int> c_and_mod_M_to_d;

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
          c_and_mod_M_to_d.insert(std::make_pair(std::make_pair(d, x), c));
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

  std::map<ManinSymbol, int> generator_to_index;
  for (auto mg: out) {
    generator_to_index.insert(std::make_pair((ManinSymbol) mg, mg.index));
  }

  return GeneratorComputationResult(out, generator_to_index, c_and_mod_M_to_d);
}

GeneratorComputationResult& compute_manin_generators(const int64_t level) {
  // [ ]: figure out a way to share this variable across translation units?
  static CacheDecorator<GeneratorComputationResult, const int64_t> _cache_compute_manin_generators(_impl_compute_manin_generators);
  return _cache_compute_manin_generators(level);
}

std::vector<ManinGenerator>& manin_generators(const int64_t level) {
  return compute_manin_generators(level).generators;
}

// Base implementation of `find_generator_index()`
// ManinGenerator _impl_find_generator(const ManinSymbol ms) {
//   auto& generators = manin_generators(ms.N);
//   auto first = generators.begin();
//   auto last = generators.end();

//   auto mg = std::find_if(first, last,
//     [&](ManinSymbol gen) { return gen.is_equivalent(ms); }
//   );

//   assert (mg != last); // Manin symbol should match one of the generators
//   return *mg;
// }

ManinGenerator find_generator(const ManinSymbol ms) {
  // static CacheDecorator<ManinGenerator, const ManinSymbol> _cache_find_generator(_impl_find_generator);
  // return _cache_find_generator(ms);
  auto& res = compute_manin_generators(ms.N);

  if (ms.c % ms.N == 0) return res.generators[0];

  fmpz_t A, B, C, D, M, N, I;

  fmpz_init_set_si(A, ms.c);
  fmpz_init_set_si(B, ms.d);
  fmpz_init_set_si(N, ms.N);

  fmpz_init(C);
  fmpz_init(D);
  fmpz_init(M);
  fmpz_init(I);

  fmpz_gcd(C, A, N);
  fmpz_divexact(A, A, C);
  fmpz_divexact(M, N, C);
  fmpz_invmod(I, A, M);
  fmpz_mul(D, B, I);
  fmpz_mod(D, D, M);

  int c = fmpz_get_si(C);
  int d_mod_M = fmpz_get_si(D);

  fmpz_clear(A);
  fmpz_clear(B);
  fmpz_clear(C);
  fmpz_clear(D);
  fmpz_clear(M);
  fmpz_clear(N);
  fmpz_clear(I);

  int d = res.c_and_mod_M_to_d[std::make_pair(c, d_mod_M)];
  ManinSymbol out_ms = {.N = ms.N, .c = c, .d = d};

  int index = res.generator_to_index[out_ms];
  return ManinGenerator(index, out_ms);
}

ManinGenerator ManinSymbol::as_generator() {
  return find_generator(*this);
}