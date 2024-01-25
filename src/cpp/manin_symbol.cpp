#include <stdio.h>
#include <vector>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fexpr_builtin.h>
#include <flint/fmpz_poly.h>
#include <flint/arith.h>

typedef struct {
  int64_t c;
  int64_t d;
  int64_t N;
} manin_generator_t;

void manin_generator_print(manin_generator_t* mg) {
  printf("(%lld, %lld)_%lld", mg->c, mg->d, mg->N);
}

typedef struct {
  manin_generator_t generator;
  fmpq_t coeff;
} mgwc_t; // manin generator with coefficient

typedef struct {
  int64_t N;
  std::vector<mgwc_t> components;
} manin_element_t;

typedef struct {
  int64_t a;
  int64_t b;
  int64_t c;
  int64_t d;
} modular_symbol_t;

typedef struct {
  int64_t a;
  int64_t b;
  int64_t c;
  int64_t d;
} sl2z_element_t;

std::vector<manin_generator_t> manin_generators(int64_t n) {
  fmpz_t N;
  fmpz_init_set_ui(N, n);

  fmpz_poly_t divisors; // acts as a list of divisors of N
  fmpz_poly_init(divisors);
  arith_divisors(divisors, N);

  int64_t tau = fmpz_poly_length(divisors);

  manin_generator_t first = {.N = n, .c = 0, .d = 1};

  std::vector<manin_generator_t> out = {first};

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
          manin_generator_t new_gen = {.N = n, .c = d, .d = c};
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

int main() {
  // manin_generator_t mg = {.N = 11, .c = 1, .d = 0};
  // manin_generator_print(&mg);
  std::vector<manin_generator_t> mgs = manin_generators(2 * 3 * 5 * 7 * 11 * 13 * 17);
  printf("%zu\n", mgs.size());
  return 0;
}