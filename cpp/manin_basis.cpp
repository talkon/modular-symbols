#include "manin_symbol.h"
#include "manin_basis.h"

#include "cache_decorator.h"

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include <iostream>
#include <vector>

void MGWC::print() {
  fmpq_print(&coeff);
  printf(" * [%lld]", index);
}

void ManinElement::print() {
  printf("level: %lld, components:", N);
  for (MGWC component : components) {
    printf(" + ");
    component.print();
  }
}

void ManinElement::print_with_generators() {
  for (MGWC component : components) {
    printf(" + ");
    fmpq_print(&component.coeff);
    printf(" * ");
    manin_generators(N)[component.index].print();
  }
}

// Auxiliary type to store the result of computing the Manin basis.
struct BasisComputationResult {

  // A vector of generators that form the Manin basis.
  // NOTE: think about whether we actually need the full generator here.
  std::vector<ManinGenerator> basis;

  // A vector containing the representation of each generator
  // as a linear combination of basis elements.
  std::vector<ManinElement> generator_to_basis;

  // A map from generator index to the index of the vector above.
  // This additional indirection is used to save memory used for
  // duplicate generators.
  std::vector<int64_t> generator_index_to_GTB_index;
};

BasisComputationResult _impl_compute_manin_basis(const int64_t level) {
  printf("[info] started computation of Manin basis for level %lld\n", level);
  std::vector<ManinGenerator> generators = manin_generators(level);
  int64_t num_generators = generators.size();
  printf("[info] finished computing generators\nnum_generators: %lld\n", num_generators);

  // Modulo out by Eta relations (see Cremona Ch 2.5)

  // Potential optimizations
  // [ ] move eta relations to manin_generators()
  // [ ] change ManinSymbol.is_equivalent to allow N | cd' + dc'

  // Manin generators modulo eta relations
  std::vector<ManinGenerator> filt_generators;

  // Current size of `filt_generators`
  int64_t filt_index = 0;

  // Mapping from index of generator to index in `filt_generators`
  // A value of -1 means mapped index has not been computed yet.
  std::vector<int64_t> generator_to_filt_generators (num_generators, -1);

  // printf("[info] eta relations:\n");
  for (ManinGenerator generator : generators) {
    if (generator_to_filt_generators[generator.index] == -1) {
      ManinGenerator generator_eta = generator.apply_eta().as_generator();
      // generator.print();
      // printf(", eta: ");
      // generator_eta.print();
      // printf("\n");

      generator_to_filt_generators[generator.index] = filt_index;
      generator_to_filt_generators[generator_eta.index] = filt_index;

      filt_generators.push_back(generator);
      filt_index++;
    }
  }

  int64_t num_filt_gens = filt_generators.size();
  printf("num_filt_gens: %lld\n", num_filt_gens);
  // for (int i = 0; i < num_filt_gens; i++) {
  //   filt_generators[i].print();
  //   printf(" %lld\n", filt_generators[i].index);
  // }

  // Compute S and T relations (see Stein Ch 3.3)

  // BUG: contains duplicate rows.
  // [ ] Start with filt_generators instead of generators to avoid redundant relations.
  // Also see Cremona Ch 2.5.

  bool done_S[num_generators];
  for (int i = 0; i < num_generators; i++) {
    done_S[i] = false;
  }
  std::vector<std::vector<int64_t>> S_rows;

  // printf("[info] S relations:\n");
  for (ManinGenerator generator : generators) {
    if (!done_S[generator.index]) {
      ManinGenerator generator_S = generator.apply_S().as_generator();
      // generator.print();
      // printf(", S: ");
      // generator_S.print();
      // printf("\n");

      std::vector<int64_t> row (num_filt_gens, 0);
      row[generator_to_filt_generators[generator.index]]++;
      row[generator_to_filt_generators[generator_S.index]]++;
      S_rows.push_back(row);

      done_S[generator.index] = true;
      done_S[generator_S.index] = true;
    }
  }

  int64_t num_S_rows = S_rows.size();

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

  printf("[info] finished computing relations\n");

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

  // fmpz_mat_print_pretty(ST);
  printf("nrows: %ld\nncols: %ld\n", fmpz_mat_nrows(ST), fmpz_mat_ncols(ST));

  fmpz_t den;
  fmpz_init(den);
  int64_t rank = fmpz_mat_rref(ST, den, ST);
  int64_t basis_size = fmpz_mat_ncols(ST) - rank;

  printf("[info] finished computing rref\n");
  printf("rref denom: ");
  fmpz_print(den);
  // printf("\nrref: \n");
  // fmpz_mat_print_pretty(ST);
  printf("\n");
  printf("rank: %lld, basis_size: %lld\n", rank, basis_size);

  // Extract information from the RREF form
  fmpz_t neg_den;
  fmpz_init(neg_den);
  fmpz_neg(neg_den, den);

  std::vector<ManinGenerator> basis;
  std::vector<ManinElement> generator_to_basis;

  int previous_pivot = -1;
  // rank + 1 is needed to add the basis elements at the end.
  for (int row = 0; row < rank + 1; row++) {
    int pivot_index;
    for (int col = previous_pivot + 1; col < num_filt_gens; col++) {
      // If nonzero element in row, this is the next pivot.
      if (!(fmpz_is_zero(fmpz_mat_entry(ST, row, col)))) {
        pivot_index = col;
        break;
      }
      // Everything else is a nonpivot (and thus in the basis)
      ManinGenerator generator = generators[filt_generators[col].index];
      basis.push_back(generator);
      generator_to_basis.push_back(generator.as_element_unchecked());
    }
    // We found a pivot, so now we construct the ManinElement corresponding to that pivot.
    previous_pivot = pivot_index;
    std::vector<MGWC> components;
    for (int col = pivot_index + 1; col < num_filt_gens; col++) {
      if (!(fmpz_is_zero(fmpz_mat_entry(ST, row, col)))) {
        fmpq_t coeff;
        fmpq_init(coeff);
        fmpq_set_fmpz_frac(coeff, fmpz_mat_entry(ST, row, col), neg_den);
        components.push_back({.index = filt_generators[col].index, .coeff = *coeff});
        fmpq_clear(coeff);
      }
    }
    ManinElement element = {.N = level, .components = components};
    generator_to_basis.push_back(element);
    // printf("[%lld] = ", filt_generators[pivot_index].index);
    // element.print();
    // printf("\n");
  }

  fmpz_mat_clear(ST);
  fmpz_clear(den);
  fmpz_clear(neg_den);

  // printf("[info] lengths\n");
  // printf(".basis: %zu\n", basis.size());
  // printf(".generator_to_basis: %zu\n", generator_to_basis.size());
  // printf(".generator_index_to_GTB_index: %zu\n", generator_to_filt_generators.size());

  return {
    .basis = basis,
    .generator_to_basis = generator_to_basis,
    .generator_index_to_GTB_index = generator_to_filt_generators
  };
}

BasisComputationResult compute_manin_basis(const int64_t level) {
  static CacheDecorator<BasisComputationResult, const int64_t> _cache_compute_manin_basis(_impl_compute_manin_basis);
  return _cache_compute_manin_basis(level);
}

ManinElement level_and_index_to_basis(const int64_t level, const int64_t index) {
  const BasisComputationResult result = compute_manin_basis(level);
  const int64_t GTB_index = result.generator_index_to_GTB_index[index];
  return result.generator_to_basis[GTB_index];
}

std::vector<ManinGenerator> manin_basis(const int64_t level) {
  const BasisComputationResult result = compute_manin_basis(level);
  return result.basis;
}