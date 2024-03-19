#include "manin_symbol.h"
#include "manin_basis.h"
#include "manin_element.h"

#include "debug_utils.h"
#include "cache_decorator.h"

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include <iostream>
#include <vector>

ManinElement ManinBasisElement::as_element() {
  fmpq_t one;
  fmpq_init(one);
  fmpq_one(one);
  std::vector<MBEWC> components = {MBEWC(basis_index, one)};
  fmpq_clear(one);

  ManinElement result = ManinElement(N, components);
  result.mark_as_sorted_unchecked();
  return result;
}

void ManinBasisElement::print_with_indices() const {
  this->print();
  printf(" [i: %lld, bi: %lld]", index, basis_index);
}

// Auxiliary type to store the result of computing the Manin basis.
struct BasisComputationResult {

  // A vector of generators that form the Manin basis.
  // NOTE: think about whether we actually need the full generator here.
  std::vector<ManinBasisElement> basis;

  // A vector containing the representation of each generator
  // as a linear combination of basis elements.
  std::vector<ManinElement> generator_to_basis;

  // A map from generator index to the index of the vector above.
  // This additional indirection is used to save memory used for
  // duplicate generators.
  std::vector<int64_t> generator_index_to_GTB_index;
};

BasisComputationResult _impl_compute_manin_basis(const int64_t level) {
  DEBUG_INFO_PRINT(2, "Started computation of Manin basis for level %lld\n", level);
  std::vector<ManinGenerator> generators = manin_generators(level);
  int64_t num_generators = generators.size();
  DEBUG_INFO_PRINT(2, "Finished computing generators, num_generators: %lld\n", num_generators);

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

  DEBUG_INFO_PRINT(5, "eta relations:\n")

  for (ManinGenerator generator : generators) {
    if (generator_to_filt_generators[generator.index] == -1) {
      ManinGenerator generator_eta = generator.apply_eta().as_generator();

      DEBUG_INFO(5,
        {
          generator.print();
          printf(", eta: ");
          generator_eta.print();
          printf("\n");
        }
      );

      generator_to_filt_generators[generator.index] = filt_index;
      generator_to_filt_generators[generator_eta.index] = filt_index;

      filt_generators.push_back(generator);
      filt_index++;
    }
  }

  int64_t num_filt_gens = filt_generators.size();
  DEBUG_INFO_PRINT(3, "num_filt_gens: %lld\n", num_filt_gens);
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

  DEBUG_INFO_PRINT(5, " S relations:\n")

  for (ManinGenerator generator : generators) {
    if (!done_S[generator.index]) {
      ManinGenerator generator_S = generator.apply_S().as_generator();

      DEBUG_INFO(5,
        {
          generator.print();
          printf(", S: ");
          generator_S.print();
          printf("\n");
        }
      )

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


  DEBUG_INFO_PRINT(3, "Finished computing relations\nnrows: %ld\nncols: %ld\n", fmpz_mat_nrows(ST), fmpz_mat_ncols(ST));

  fmpz_t den;
  fmpz_init(den);
  int64_t rank = fmpz_mat_rref(ST, den, ST);
  int64_t basis_size = fmpz_mat_ncols(ST) - rank;

  DEBUG_INFO(3,
    {
      printf("Finished computing rref\n");
      printf("rref denom: ");
      fmpz_print(den);
      // printf("\nrref: \n");
      // fmpz_mat_print_pretty(ST);
      printf("\n");
      printf("rank: %lld, basis_size: %lld\n", rank, basis_size);
    }
  )

  // Extract information from the RREF form
  fmpz_t neg_den;
  fmpz_init(neg_den);
  fmpz_neg(neg_den, den);

  std::map<int, ManinBasisElement> basis_map;
  std::vector<ManinElement> generator_to_basis;

  std::map<int64_t, int64_t> index_to_basis_index;
  int64_t bi = 0;

  auto get_basis_index = [&index_to_basis_index, &bi] (int64_t i) {
    if (auto it = index_to_basis_index.find(i); it != index_to_basis_index.end()) {
      return it->second;
    } else {
      index_to_basis_index.insert(std::make_pair(i, bi));
      bi++;
      return bi - 1;
    }
  };

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
      int64_t basis_index = get_basis_index(filt_generators[col].index);
      ManinBasisElement mbe(basis_index, generator);
      basis_map.insert(std::pair(basis_index, mbe));
      generator_to_basis.push_back(mbe.as_element());
    }
    // We found a pivot, so now we construct the ManinElement corresponding to that pivot.
    previous_pivot = pivot_index;
    std::vector<MBEWC> components;
    for (int col = pivot_index + 1; col < num_filt_gens; col++) {
      if (!(fmpz_is_zero(fmpz_mat_entry(ST, row, col)))) {
        fmpq_t coeff;
        fmpq_init(coeff);
        fmpq_set_fmpz_frac(coeff, fmpz_mat_entry(ST, row, col), neg_den);
        int64_t basis_index = get_basis_index(filt_generators[col].index);
        components.push_back(MBEWC(basis_index, coeff));
        fmpq_clear(coeff);
      }
    }
    ManinElement element = ManinElement(level, components);
    element.sort();
    generator_to_basis.push_back(element);
    // printf("[%lld] = ", filt_generators[pivot_index].index);
    // element.print();
    // printf("\n");
  }

  std::vector<ManinBasisElement> basis;
  for (int bi = 0; bi < basis_size; bi++) {
    auto pair = basis_map.find(bi);
    basis.push_back(pair->second);
  }

  fmpz_mat_clear(ST);
  fmpz_clear(den);
  fmpz_clear(neg_den);

  DEBUG_INFO(5,
    {
      printf(" basis computation result\n");
      printf(".basis: %zu\n", basis.size());
      for (auto mbe: basis) {
        mbe.print_with_indices();
        printf("\n");
      }
      printf(".generator_to_basis: %zu\n", generator_to_basis.size());
      for (int i = 0; i < generator_to_basis.size(); i++) {
        generators[i].print();
        printf(" = ");
        generator_to_basis[i].print();
        printf("\n");
      }
      printf(".generator_index_to_GTB_index: %zu\n", generator_to_filt_generators.size());
      for (auto i: generator_to_filt_generators) {
        printf("%lld ", i);
      }
      printf("\n");
    }
  )

  DEBUG_INFO_PRINT(2, "Finished computation of Manin basis for level %lld\n", level);

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

std::vector<ManinBasisElement> manin_basis(const int64_t level) {
  const BasisComputationResult result = compute_manin_basis(level);
  return result.basis;
}

ManinElement fraction_to_manin_element(const int64_t a, const int64_t b, const int64_t level) {

  // printf("fraction_to_manin_element called with (a, b, level) = (%lld, %lld, %lld)\n", a, b, level);

  if (a == 0) {
    return ManinElement::zero(level);
  }

  if (b == 0) {
    // (1, 0)_N is always the first element in the basis.
    // XXX: actually check this
    ManinSymbol ms = {.c = 1, .d = 0, .N = level};
    ManinGenerator mg = find_generator(ms);
    ManinBasisElement mbe(0, mg);
    return mbe.as_element();
  }

  ManinElement result = ManinElement::zero(level);

  // XXX: Usage of fmpz might not be necessary -- I think everything should fit in 64 bits.
  fmpz_t X, Y, Xp, Xq, Yp, Yq, Q, R;
  fmpz_init(Q);
  fmpz_init(R);

  fmpz_init_set_si(X, a);
  fmpz_init_set_si(Y, b);

  fmpz_init_set_si(Xp, 0);
  fmpz_init_set_si(Xq, 1);
  fmpz_init_set_si(Yp, 1);
  fmpz_init_set_si(Yq, 0);

  int iter = 0;

  while (true) {
    // X = Y * Q + R, or R = X - Y * Q
    fmpz_tdiv_qr(Q, R, X, Y);
    fmpz_set(X, R);
    fmpz_submul(Xp, Yp, Q);
    fmpz_submul(Xq, Yq, Q);
    iter++;

    if (iter > 1) {
      int64_t xq = fmpz_get_si(Xq);
      int64_t yq = fmpz_get_si(Yq);
      ManinSymbol ms = {.N = level, .c = -xq, .d = -yq};
      ManinElement me = level_and_index_to_basis(level, ms.as_generator().index);

      // ms.print();
      // printf(" = ");
      // me.print_with_generators();
      // printf("\n");

      result -= me;
    }

    if (fmpz_is_zero(X)) {
      break;
    }

    // Y = X * Q + R, or R = Y - X * Q
    fmpz_tdiv_qr(Q, R, Y, X);
    fmpz_set(Y, R);
    fmpz_submul(Yp, Xp, Q);
    fmpz_submul(Yq, Xq, Q);
    iter++;

    int64_t xq = fmpz_get_si(Xq);
    int64_t yq = fmpz_get_si(Yq);
    ManinSymbol ms = {.N = level, .c = -yq, .d = xq};
    ManinElement me = level_and_index_to_basis(level, ms.as_generator().index);

    // ms.print();
    // printf(" = ");
    // me.print_with_generators();
    // printf("\n");

    result -= me;

    if (fmpz_is_zero(Y)) {
      break;
    }
  }

  fmpz_clear(X);
  fmpz_clear(Y);
  fmpz_clear(Xp);
  fmpz_clear(Xq);
  fmpz_clear(Yp);
  fmpz_clear(Yq);
  fmpz_clear(Q);
  fmpz_clear(R);

  return result;
}