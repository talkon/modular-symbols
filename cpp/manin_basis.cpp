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
#include <map>
#include <queue>
#include <set>
#include <cassert>

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

  BasisComputationResult(
    std::vector<ManinBasisElement> basis,
    std::vector<ManinElement> generator_to_basis,
    std::vector<int64_t> generator_index_to_GTB_index
  ) :
    basis(basis),
    generator_to_basis(generator_to_basis),
    generator_index_to_GTB_index(generator_index_to_GTB_index)
  {}
};

BasisComputationResult _impl_compute_manin_basis(const int64_t level) {
  DEBUG_INFO_PRINT(3, "Started computation of Manin basis for level %lld\n", level);
  std::vector<ManinGenerator> generators = manin_generators(level);
  int64_t num_generators = generators.size();
  DEBUG_INFO_PRINT(3, "Finished computing generators, num_generators: %lld\n", num_generators);

  // ------------------------------------------------ //
  // Modulo out by Eta relations (see Cremona Ch 2.5) //
  // ------------------------------------------------ //

  // Potential optimizations
  // [ ] move eta relations to manin_generators()
  // [ ] change ManinSymbol.is_equivalent to allow N | cd' + dc'

  // Manin generators modulo eta relations
  std::vector<ManinGenerator> eta_generators;

  // Current size of `eta_generators`
  int64_t eta_index = 0;

  // Mapping from index of generator to index in `eta_generators`
  // A value of -1 means mapped index has not been computed yet.
  std::vector<int64_t> generator_to_eta_generators (num_generators, -1);

  DEBUG_INFO_PRINT(6, "eta relations:\n")

  for (ManinGenerator generator : generators) {
    if (generator_to_eta_generators[generator.index] == -1) {
      ManinGenerator generator_eta = generator.apply_eta().as_generator();

      DEBUG_INFO(7,
        {
          generator.print();
          printf(", eta: ");
          generator_eta.print();
          printf("\n");
        }
      );

      generator_to_eta_generators[generator.index] = eta_index;
      generator_to_eta_generators[generator_eta.index] = eta_index;

      eta_generators.push_back(generator);
      eta_index++;
    }
  }

  int64_t num_eta_gens = eta_generators.size();
  assert(eta_index = num_eta_gens);

  DEBUG_INFO_PRINT(4, "num_eta_gens: %lld\n", num_eta_gens);

  DEBUG_INFO(7,
    {
      for (int i = 0; i < num_eta_gens; i++) {
        printf("%d ", i);
        eta_generators[i].print();
        printf(" %lld\n", eta_generators[i].index);
      }
    }
  )

  // -------------------------------------------- //
  // Modulo out by S relations (see Stein Ch 3.3) //
  // -------------------------------------------- //

  // BUG: contains duplicate rows.
  // [ ] Start with eta_generators instead of generators to avoid redundant relations.
  // Also see Cremona Ch 2.5.

  bool done_S[num_eta_gens];
  for (int i = 0; i < num_eta_gens; i++) {
    done_S[i] = false;
  }

  // Manin generators modulo eta and S generators
  std::vector<ManinGenerator> S_generators_pos;
  std::vector<ManinGenerator> S_generators_neg;

  // Current index of `S_generators`
  int64_t S_index = 1;

  // Mapping from index of generator to index in `eta_generators`
  // A value of 0 means mapped index has not been computed yet, or
  // the eta generator is zero
  std::vector<int64_t> eta_to_S_generators (num_eta_gens, 0);

  for (ManinGenerator eta_gen : eta_generators) {

    int index = generator_to_eta_generators[eta_gen.index];

    if (!done_S[index]) {

      ManinGenerator eta_gen_S = eta_gen.apply_S().as_generator();
      int index_S = generator_to_eta_generators[eta_gen_S.index];
      eta_gen_S = eta_generators[index_S];

      if (eta_gen_S.index == eta_gen.index) {

        // eta_gen.print();
        // printf(" -> zero\n");

        done_S[index] = true;
        eta_to_S_generators[eta_gen.index] = 0;

      } else {

        done_S[index] = true;
        done_S[index_S] = true;

        eta_to_S_generators[index] = S_index;
        eta_to_S_generators[index_S] = -S_index;

        S_generators_pos.push_back(eta_gen);
        S_generators_neg.push_back(eta_gen_S);

        S_index++;

      }
    }
  }

  int64_t num_S_gens = S_generators_pos.size();
  assert(S_index - 1 == num_S_gens);
  assert(S_index - 1 == S_generators_neg.size());


  std::vector<int64_t> generator_to_S_generators (num_generators, 0);
  for (int i = 0; i < num_generators; i++) {
    generator_to_S_generators[i] = eta_to_S_generators[generator_to_eta_generators[i]];
  }

  DEBUG_INFO_PRINT(4, "num_S_gens: 2 * %lld\n", num_S_gens);

  DEBUG_INFO(7,
    {
      for (int i = 0; i < num_S_gens; i++) {
        S_generators_pos[i].print();
        int idx = S_generators_pos[i].index;
        printf(" id %d eta %lld s %lld \n", idx, generator_to_eta_generators[idx], generator_to_S_generators[idx]);
      }

      for (int i = 0; i < num_S_gens; i++) {
        S_generators_neg[i].print();
        int idx = S_generators_neg[i].index;
        printf(" id %d eta %lld s %lld \n", idx, generator_to_eta_generators[idx], generator_to_S_generators[idx]);
      }
    }
  )

  // -------------------------------------------- //
  // Compute T relation matrix (see Stein Ch 3.3) //
  // -------------------------------------------- //

  // std::vector<ManinGenerator> S_generators = S_generators_pos;
  // S_generators.insert(S_generators.end(), S_generators_neg.begin(), S_generators_neg.end());

  // bool done_T_pos[num_S_gens];
  // bool done_T_neg[num_S_gens];
  // for (int i = 0; i < num_S_gens; i++) {
  //   done_T_pos[i] = false;
  //   done_T_neg[i] = false;
  // }

  bool done_T[num_eta_gens];
  for (int i = 0; i < num_eta_gens; i++) {
    done_T[i] = false;
  }

  std::set<std::tuple<int64_t, int64_t, int64_t>> T_orbits;
  std::vector<std::vector<int64_t>> T_rows;

  // XXX: think about this more and figure out what the math actually is.
  for (ManinGenerator generator : generators) {

    int index = generator_to_S_generators[generator.index];
    int eta_index = generator_to_eta_generators[generator.index];

    // if (index == 0) {
    //   throw std::runtime_error("S_generators should not contain zero elements");
    // }

    // bool done = done_T[eta_index];
    // // bool done = false;

    // if (!done) {

    ManinGenerator generator_T = generator.apply_T().as_generator();
    ManinGenerator generator_T_2 = generator.apply_T_2().as_generator();

    int index_T = generator_to_S_generators[generator_T.index];
    int index_T_2 = generator_to_S_generators[generator_T_2.index];

    DEBUG_INFO(7,
      {
        printf("%4lld %4lld %4lld -> ", generator.index, generator_T.index, generator_T_2.index);
        printf("%4lld %4lld %4lld -> ",
          generator_to_eta_generators[generator.index],
          generator_to_eta_generators[generator_T.index],
          generator_to_eta_generators[generator_T_2.index]
        );
        printf("%4d %4d %4d\n",
          index,
          index_T,
          index_T_2
        );
      }
    )

    if ((index == 0 || index_T == 0 || index_T_2 == 0) && (index + index_T + index_T_2 == 0)) continue;

    std::tuple<int64_t, int64_t, int64_t> orbit = {0, 0, 0};

    if (index <= index_T && index_T <= index_T_2)   orbit = {index, index_T, index_T_2};
    if (index <= index_T_2 && index_T_2 <= index_T) orbit = {index, index_T_2, index_T};
    if (index_T <= index && index <= index_T_2)     orbit = {index_T, index, index_T_2};
    if (index_T <= index_T_2 && index_T_2 <= index) orbit = {index_T, index_T_2, index};
    if (index_T_2 <= index && index <= index_T)     orbit = {index_T_2, index, index_T};
    if (index_T_2 <= index_T && index_T <= index)   orbit = {index_T_2, index_T, index};

    if (orbit == std::make_tuple(0, 0, 0)) throw std::runtime_error("orbit should not be trivial");

    std::tuple<int64_t, int64_t, int64_t> negated_orbit = {
      -std::get<2>(orbit),
      -std::get<1>(orbit),
      -std::get<0>(orbit)
    };

    DEBUG_INFO_PRINT(6, "orbit: %lld %lld %lld, negated_orbit: %lld %lld %lld\n",
      std::get<0>(orbit),
      std::get<1>(orbit),
      std::get<2>(orbit),
      std::get<0>(negated_orbit),
      std::get<1>(negated_orbit),
      std::get<2>(negated_orbit)
    );

    if (auto search = T_orbits.find(orbit); search != T_orbits.end()) continue;
    if (auto search = T_orbits.find(negated_orbit); search != T_orbits.end()) continue;

    T_orbits.insert(orbit);

    std::vector<int64_t> row (num_S_gens, 0);

    // bool mark_as_done = !(index == 0 || index_T == 0 || index_T_2 == 0);

    // done_T[generator_to_eta_generators[generator.index]] |= mark_as_done;
    // done_T[generator_to_eta_generators[generator_T.index]] |= mark_as_done;
    // done_T[generator_to_eta_generators[generator_T_2.index]] |= mark_as_done;

    for (int idx : {index, index_T, index_T_2}) {

      if (idx > 0) {
        row[idx - 1]++;
        // done_T_pos[idx - 1] |= mark_as_done;
      }

      if (idx < 0) {
        row[-idx - 1]--;
        // done_T_neg[-idx - 1] |= mark_as_done;
      }

    }

    DEBUG_INFO(7,
      {
        printf("new orbit!\n");
        printf("row: ");
        for (auto x : row) printf("%lld ", x);
        printf("\n");
      }
    )

    T_rows.push_back(row);
  }

  int64_t num_T_rows = T_rows.size();

  fmpz_mat_t T_mat;
  fmpz_mat_init(T_mat, num_T_rows, num_S_gens);

  for (int row = 0; row < num_T_rows; row++) {
    for (int col = 0; col < num_S_gens; col++) {
      fmpz_set_si(fmpz_mat_entry(T_mat, row, col), T_rows[row][col]);
    }
  }

  DEBUG_INFO_PRINT(4, "Finished computing relations\nnrows: %ld\nncols: %ld\n", fmpz_mat_nrows(T_mat), fmpz_mat_ncols(T_mat));


  DEBUG_INFO(7,
    {
      printf("T_mat: \n");
      fmpz_mat_print_pretty(T_mat);
    }
  )

#define NEW_RR 1
#if NEW_RR

  // ------------------------------ //
  // Row reducing T relation matrix //
  // ------------------------------ //

  // 0. Preprocess T relation matrix
  // If we see any +-3, this corresponds to a fixed point of T (i.e. will be zero), so we can replace it by 1
  // TODO: this is a bit hacky, think about what's actually happening here and make sure we're not missing any other case (like a col with 1 and 2).
  for (int row = 0; row < num_T_rows; row++) {
    for (int col = 0; col < num_S_gens; col++) {
      if (fmpz_equal_si(fmpz_mat_entry(T_mat, row, col), 3)) {
        fmpz_set_si(fmpz_mat_entry(T_mat, row, col), 1);
      }
      if (fmpz_equal_si(fmpz_mat_entry(T_mat, row, col), -3)) {
        fmpz_set_si(fmpz_mat_entry(T_mat, row, col), 1);
      }
    }
  }

  // 1. Constructing graph G with rows as vertices and columns as edges

  // r_1 -> r_2 -> col
  std::vector<std::map<int, int>> G(num_T_rows);

  // (col, row) such that row has the only nonzero entry in col
  std::vector<std::pair<int, int>> single_cols;

  for (int col = 0; col < num_S_gens; col++) {
    int nonzero_1 = -1;
    int nonzero_2 = -1;
    for (int row = 0; row < num_T_rows; row++) {
      if (!fmpz_is_zero(fmpz_mat_entry(T_mat, row, col))) {
        if (nonzero_1 == -1) nonzero_1 = row;
        else nonzero_2 = row;
      }
    }

    if (nonzero_1 == -1) continue;
    else if (nonzero_2 == -1) single_cols.emplace_back(col, nonzero_1);
    else {
      G[nonzero_1].insert(std::make_pair(nonzero_2, col));
      G[nonzero_2].insert(std::make_pair(nonzero_1, col));
    }
  }

  // 2. Compute BFS tree
  assert (single_cols.size() > 0);
  int start_row = single_cols[0].second;
  int start_col = single_cols[0].first;

  std::vector<bool> seen(num_T_rows, false);
  std::queue<int> worklist;

  // stores edges row1 -> row2 in tree as (col, row1, row2)
  std::vector<std::tuple<int, int, int>> tree;

  worklist.push(start_row);
  seen[start_row] = true;

  while (worklist.size() > 0) {

    int node = worklist.front();
    worklist.pop();

    auto adj = G[node];
    for (auto const& [r, c] : adj) {
      if (!seen[r]) {
        tree.emplace_back(c, node, r);
        worklist.push(r);
        seen[r] = true;
      }
    }
  }

  // 3. Apply row operations
  std::vector<bool> is_pivot(num_S_gens, false);
  is_pivot[start_col] = true;

  assert(tree.size() == num_T_rows - 1);
  for (auto it = tree.rbegin(); it != tree.rend(); it++) {
    auto [c, r1, r2] = *it;
    // printf("%d %d %d\n", c, r1, r2);
    is_pivot[c] = true;

    // If there are two nonzero elements in a col, then they must be -1 or +1
    int sub = fmpz_equal(fmpz_mat_entry(T_mat, r1, c), fmpz_mat_entry(T_mat, r2, c));

    if (sub) {
      // row1 <- row1 - row2
      for (int col = 0; col < num_S_gens; col++) {
        fmpz_sub(fmpz_mat_entry(T_mat, r1, col), fmpz_mat_entry(T_mat, r1, col), fmpz_mat_entry(T_mat, r2, col));
      }
    } else {
      // row1 <- row1 + row2
      for (int col = 0; col < num_S_gens; col++) {
        fmpz_add(fmpz_mat_entry(T_mat, r1, col), fmpz_mat_entry(T_mat, r1, col), fmpz_mat_entry(T_mat, r2, col));
      }
    }
  }

  DEBUG_INFO(7,
    {
      printf("T_mat, after RR: \n");
      fmpz_mat_print_pretty(T_mat);
    }
  )

  // 4. Construct output
  std::vector<ManinBasisElement> basis;
  int basis_index = 0;
  std::map<int, int> col_to_basis_index;
  std::map<int, ManinElement> generator_to_basis;

  // Make basis (i.e. nonpivots)
  for (int col = 0; col < num_S_gens; col++) {
    if (!is_pivot[col]) {
      ManinGenerator generator = S_generators_pos.at(col);
      ManinBasisElement mbe(basis_index, generator);
      basis.push_back(mbe);
      col_to_basis_index.insert(std::make_pair(col, basis_index));
      generator_to_basis.insert(std::make_pair(col, mbe.as_element()));
      basis_index++;
    }
  }

  // Make GTB for pivots
  tree.emplace_back(start_col, -1, start_row);
  for (auto it = tree.begin(); it != tree.end(); it++) {
    auto [c, r1, r2] = *it;
    // Row r2 contains pivot at column c, which is either +1 or -1

    // If pivot is +1, then we have to negate the rest
    int negate = fmpz_is_one(fmpz_mat_entry(T_mat, r2, c));

    std::vector<MBEWC> components;
    for (int col = 0; col < num_S_gens; col++) {
      if (col != c && !(fmpz_is_zero(fmpz_mat_entry(T_mat, r2, col)))) {
        fmpq_t coeff;
        fmpq_init(coeff);
        fmpq_set_fmpz(coeff, fmpz_mat_entry(T_mat, r2, col));

        if (negate) fmpq_neg(coeff, coeff);

        int64_t basis_index = col_to_basis_index[col];
        components.emplace_back(basis_index, coeff);
        fmpq_clear(coeff);
      }
    }

    ManinElement element = ManinElement(level, components);
    element.mark_as_sorted_unchecked();
    generator_to_basis.insert(std::make_pair(c, element));
  }

// ---- (old code:) ---- //
#else

  fmpz_t den;
  fmpz_init(den);

  int64_t rank = fmpz_mat_rref_mul(T_mat, den, T_mat);
  int64_t basis_size = fmpz_mat_ncols(T_mat) - rank;

  DEBUG_INFO(3,
    {
      printf("Finished computing rref\n");
      printf("rref denom: ");
      fmpz_print(den);
      printf("\nrank: %lld, basis_size: %lld\n", rank, basis_size);
    }
  )

  DEBUG_INFO(7,
    {
      printf("rref: \n");
      fmpz_mat_print_pretty(T_mat);
      printf("\n");
    }
  )

  // Extract information from the RREF form
  fmpz_t neg_den;
  fmpz_init(neg_den);
  fmpz_neg(neg_den, den);

  std::vector<ManinBasisElement> basis;
  std::map<int, int> col_to_basis_index;
  std::map<int, ManinElement> generator_to_basis;

  int basis_index = 0;

  // Find nonpivots
  int previous_pivot = -1;
  for (int row = 0; row < rank + 1; row++) {
    for (int col = previous_pivot + 1; col < num_S_gens; col++) {
      // If nonzero element in row, this is the next pivot.
      if (row < rank && !(fmpz_is_zero(fmpz_mat_entry(T_mat, row, col)))) {
        previous_pivot = col;
        break;
      }

      // Everything else is a nonpivot (and thus in the basis)
      // printf("(%d, %d)\n", row, col);
      ManinGenerator generator = S_generators_pos.at(col);

      // generator.print();
      // printf("\n");

      ManinBasisElement mbe(basis_index, generator);
      basis.push_back(mbe);
      col_to_basis_index.insert(std::make_pair(col, basis_index));
      generator_to_basis.insert(std::make_pair(col, mbe.as_element()));

      // printf("nonpivot %d -> ", col);
      // mbe.as_element().print();
      // printf("\n");

      basis_index++;
    }
  }

  // Fix: order of things added to generator_to_basis;
  // Deals with pivots
  previous_pivot = -1;
  for (int row = 0; row < rank; row++) {
    int pivot_index;
    for (int col = previous_pivot + 1; col < num_S_gens; col++) {
      // If nonzero element in row, this is the next pivot.
      if (!(fmpz_is_zero(fmpz_mat_entry(T_mat, row, col)))) {
        pivot_index = col;
        break;
      }
      // Everything else is a nonpivot (and thus in the basis)
    }
    // We found a pivot, so now we construct the ManinElement corresponding to that pivot.
    previous_pivot = pivot_index;
    std::vector<MBEWC> components;
    for (int col = pivot_index + 1; col < num_S_gens; col++) {
      if (!(fmpz_is_zero(fmpz_mat_entry(T_mat, row, col)))) {
        fmpq_t coeff;
        fmpq_init(coeff);
        fmpq_set_fmpz_frac(coeff, fmpz_mat_entry(T_mat, row, col), neg_den);
        int64_t basis_index = col_to_basis_index[col];
        components.emplace_back(basis_index, coeff);
        fmpq_clear(coeff);
      }
    }
    ManinElement element = ManinElement(level, components);
    element.mark_as_sorted_unchecked();
    generator_to_basis.insert(std::make_pair(pivot_index, element));

    // printf("pivot %d -> ", pivot_index);
    // element.print();
    // printf("\n");

  }

  fmpz_clear(den);
  fmpz_clear(neg_den);

#endif

  fmpz_mat_clear(T_mat);

  std::vector<ManinElement> gtb_vec;
  for (int i = 0; i < num_S_gens; i++) {
    gtb_vec.push_back(generator_to_basis.at(i));
  }

  DEBUG_INFO(6,
    {
      printf(" basis computation result\n");
      printf(".basis: %zu\n", basis.size());
      for (auto mbe: basis) {
        mbe.print_with_indices();
        printf("\n");
      }
      printf(".generator_to_basis: %zu\n", generator_to_basis.size());
      for (int i = 0; i < gtb_vec.size(); i++) {
        S_generators_pos.at(i).print();
        printf(" = ");
        gtb_vec[i].print();
        printf("\n");
      }
      printf(".generator_index_to_GTB_index: %zu\n", generator_to_S_generators.size());
      for (auto i: generator_to_S_generators) {
        printf("%lld ", i);
      }
      printf("\n");
    }
  )

  DEBUG_INFO_PRINT(3, "Finished computation of Manin basis for level %lld\n", level);

  return BasisComputationResult(basis, gtb_vec, generator_to_S_generators);
}

// FIXME: figure out a way to return cached things as references instead of values in order to avoid unnecessary copying.
BasisComputationResult& compute_manin_basis(const int64_t level) {
  static CacheDecorator<BasisComputationResult, const int64_t> _cache_compute_manin_basis(_impl_compute_manin_basis);
  return _cache_compute_manin_basis(level);
}

ManinElement level_and_index_to_basis(const int64_t level, const int64_t index) {
  const BasisComputationResult& result = compute_manin_basis(level);
  const int64_t GTB_index = result.generator_index_to_GTB_index[index];
  if (GTB_index > 0) {
    return result.generator_to_basis[GTB_index - 1];
  } else if (GTB_index == 0) {
    return ManinElement::zero(level);
  } else {
    return result.generator_to_basis[-GTB_index - 1].negate();
  }
}

std::vector<ManinBasisElement> manin_basis(const int64_t level) {
  const BasisComputationResult& result = compute_manin_basis(level);
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
    ManinSymbol ms = {.N = level, .c = 1, .d = 0};
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
