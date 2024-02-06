#include "manin_symbol.h"

#include "cache_decorator.h"

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>
#include <flint/fexpr_builtin.h>
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


// Base implementation of `manin_generators()`
// XXX: Surely there is a better approach to this pattern, but I think this works for now.
// [ ] Check that caching actually works across multiple compilation units.
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

void manin_basis(int64_t n) {
  printf("here 000\n");
  std::vector<ManinGenerator> generators = manin_generators(n);
  int64_t num_generators = generators.size();
  printf("num_generators: %lld\n", num_generators);

  printf("here 001\n");

  // Modulo out by Eta relations (see Cremona Ch 2.5)

  // Potential optimizations
  // [ ] move eta relations to manin_generators()
  // [ ] change ManinSymbol.is_equivalent to allow N | cd' + dc'

  // Manin generators modulo eta relations
  std::vector<ManinGenerator*> filt_generators;

  // Current size of `filt_generators`
  int64_t filt_index = 0;

  // Mapping from index of generator to index in `filt_generators`
  // A value of -1 means mapped index has not been computed yet.
  int64_t generator_to_filt_generators[num_generators];
  for (int i = 0; i < num_generators; i++) {
    generator_to_filt_generators[i] = -1;
  }

  printf("here 001b\n");

  for (ManinGenerator generator : generators) {
    generator.print();
    printf("\n");
    if (generator_to_filt_generators[generator.index] == -1) {
      ManinGenerator generator_eta = generator.apply_eta().as_generator();

      printf("generator: ");
      generator.print();
      printf("\n generator.apply_eta(): ");
      generator.apply_eta().print();
      printf("\n generator_eta: ");
      _impl_find_generator(generator.apply_eta()).print();
      printf("\n");
      // ManinSymbol sym = {.N = 11, .c = 1, .d = 0};
      // ManinGenerator* found = find_generator(sym);
      // sym.print(); printf(" "); found->print();
      // printf("\n");
      // sym = {.N = 11, .c = 1, .d = 1};
      // found = find_generator(sym);
      // sym.print(); printf(" "); (*found).print();
      // printf("\n");
      // sym = {.N = 11, .c = -1, .d = 1};
      // found = find_generator(sym);
      // sym.print(); printf(" "); found->print();
      // printf("\n");

      generator_to_filt_generators[generator.index] = filt_index;
      generator_to_filt_generators[generator_eta.index] = filt_index;

      filt_generators.push_back(&generator);
      filt_index++;
    }
  }

  int64_t num_filt_gens = filt_generators.size();
  printf("num_filt_gens: %lld\n", num_filt_gens);

  printf("here 002\n");

  // Compute S and T relations (see Stein Ch 3.3)

  // BUG: contains duplicate rows.
  // [ ] Start with filt_generators instead of generators to avoid redundant relations.
  // Also see Cremona Ch 2.5.

  bool done_S[num_generators];
  for (int i = 0; i < num_generators; i++) {
    done_S[i] = false;
  }
  std::vector<std::vector<int64_t>> S_rows;

  for (ManinGenerator generator : generators) {
    if (!done_S[generator.index]) {
      ManinGenerator generator_S = generator.apply_S().as_generator();

      std::vector<int64_t> row (num_filt_gens, 0);
      printf("indices: %lld %lld\n", generator.index, generator_S.index);
      row[generator_to_filt_generators[generator.index]]++;
      row[generator_to_filt_generators[generator_S.index]]++;
      S_rows.push_back(row);

      done_S[generator.index] = true;
      done_S[generator_S.index] = true;
    }
  }

  int64_t num_S_rows = S_rows.size();

  printf("here 003\n");

  bool done_T[num_generators];
  for (int i = 0; i < num_generators; i++) {
    done_T[i] = false;
  }
  std::vector<std::vector<int64_t>> T_rows;

  for (ManinGenerator generator : generators) {
    if (!done_T[generator.index]) {
      ManinGenerator generator_T = generator.apply_T().as_generator();
      ManinGenerator generator_T_2 = generator.apply_T_2().as_generator();

      std::vector<int64_t> row (num_generators, 0);
      row[generator_to_filt_generators[generator.index]]++;
      row[generator_to_filt_generators[generator_T.index]]++;
      row[generator_to_filt_generators[generator_T_2.index]]++;
      T_rows.push_back(row);

      done_T[generator.index] = true;
      done_T[generator_T.index] = true;
      done_T[generator_T_2.index] = true;
    }
  }

  int64_t num_T_rows = T_rows.size();

  printf("here 004\n");

  // Create ST relation matrix
  fmpz_mat_t ST;
  fmpz_mat_init(ST, num_S_rows + num_T_rows, num_filt_gens);

  for (int row = 0; row < num_S_rows; row++) {
    for (int col = 0; col < num_filt_gens; col++) {
      fmpz_set_ui(fmpz_mat_entry(ST, row, col), S_rows[row][col]);
    }
  }

  for (int row = 0; row < num_T_rows; row++) {
    for (int col = 0; col < num_filt_gens; col++) {
      fmpz_set_ui(fmpz_mat_entry(ST, row + num_S_rows, col), T_rows[row][col]);
    }
  }

  printf("here 005\n");

  fmpz_mat_print_pretty(ST);
  printf("nrows: %ld\nncols: %ld", fmpz_mat_nrows(ST), fmpz_mat_ncols(ST));
  fmpz_mat_clear(ST);
}

int main() {
  // // Tests manin_generators
  // std::vector<ManinGenerator> mgs;
  // for (int i = 0; i < 10; i++) {
  //   mgs = manin_generators(2 * 3 * 5 * 7 * 11 * 13);
  //   printf("%zu\n", mgs.size());
  // }
  // mgs[2 * 3 * 5 * 7 * 11 * 13 + 50000].print();
  // printf("\n");

  // // Tests find_generator
  // for (int i = 11; i < 20; i++) {
  //   ManinSymbol ms1 = {.N = i, .c = 4, .d = 7};
  //   ManinGenerator mg1 = find_generator(ms1);
  //   ms1.print();
  //   mg1.print();
  //   printf(" %lld %p\n", mg1.index, &mg1);
  // }

  // Tests relation matrix
  manin_basis(11);

  return 0;
}