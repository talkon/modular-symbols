#include "boundary_map.h"
#include "manin_element.h"
#include "utils.h"
#include "fmpz_mat_helpers.h"

#include "debug_utils.h"
#include "cache_decorator.h"

#include <flint/ulong_extras.h>
#include <flint/fmpz_mat.h>

#include <vector>
#include <cassert>

void Cusp::print() {
  printf("(%lld, %lld)_%lld", c, d, N);
}

utils::XgcdResult Cusp::xgcd() const {
  return utils::xgcd(c, d);
}

bool Cusp::operator<(const Cusp& other) const {
  return (std::tuple<int64_t,int64_t,int64_t>) {N, c, d} < (std::tuple<int64_t,int64_t,int64_t>) {other.N, other.c, other.d};
}

Cusp Cusp::negated() const {
  return {.c = -c, .d = d, .N = N};
}

bool Cusp::is_equivalent(const Cusp& other) {
  if (N != other.N)
    return false;

  // XXX: consider precomputing this?
  utils::XgcdResult this_result = this->xgcd();
  utils::XgcdResult other_result = other.xgcd();

  assert (this_result.gcd == 1);
  assert (other_result.gcd == 1);

  int64_t g = n_gcd(this->d * other.d, N);

  return (
    (this_result.a * other.d - this->d * other_result.a) % g == 0
    || (this_result.a * other.d + this->d * other_result.a) % g == 0
  );
}

// Computes the result of the boundary map on a given Manin generator.
// Returns a pair of cusps {c1, c2}; this should be treated as + c1 - c2.
std::pair<Cusp, Cusp> boundary_map(const ManinGenerator mg) {
  utils::XgcdResult result = utils::xgcd(mg.c, mg.d);
  assert (result.gcd == 1);

  Cusp pos_cusp = {.c = result.b, .d = mg.c, .N = mg.N};
  Cusp neg_cusp = {.c = -result.a, .d = mg.d, .N = mg.N};

  return std::make_pair(pos_cusp, neg_cusp);
}

DenseBasis cuspidal_manin_basis(int64_t level) {
  std::vector<ManinBasisElement> full_basis = manin_basis(level);
  int b_size = full_basis.size();
  DEBUG_INFO_PRINT(1, "Started computation of cuspidal Manin basis for level %lld\n", level);

  std::vector<Cusp> representatives;
  // XXX: not sure result_cache is any helpful.
  std::map<Cusp, int> result_cache;

  // given a cusp, returns an index of its representative in representatives
  auto find_representative_index = [&](Cusp& cusp) {

    // In cache
    if (auto it = result_cache.find(cusp); it != result_cache.end())
      return it->second;

    // Matches known representative
    int i = 0;
    for (; i < representatives.size(); i++) {
      Cusp representative = representatives[i];
      if (representative.is_equivalent(cusp)) {
        result_cache.insert(std::make_pair(cusp, i));
        return i;
      }
    }

    // New representative class
    // printf("new representative: ");
    // cusp.print();
    // printf("\n");
    representatives.push_back(cusp);
    result_cache.insert(std::make_pair(cusp, i));
    return i;
  };

  // Compute result of boundary_map, and fill in representatives
  std::vector<std::pair<int, int>> mapped_basis;
  for (ManinGenerator mg : full_basis) {
    auto pair = boundary_map(mg);
    int first = find_representative_index(pair.first);
    int second = find_representative_index(pair.second);
    mapped_basis.push_back(std::make_pair(first, second));
  }

  DEBUG_INFO_PRINT(3, "Finished computing representatives\
  \nfull_basis size: %zu\
  \nnum_representatives: %zu\n", b_size, representatives.size());

  // Construct boundary map matrix as a dense fmpz_mat_t
  fmpz_mat_t boundary_map_matrix;
  fmpz_mat_init(boundary_map_matrix, representatives.size(), b_size);
  fmpz_mat_zero(boundary_map_matrix);

  for (int gen_index = 0; gen_index < mapped_basis.size(); gen_index++) {
    int pos_index = mapped_basis[gen_index].first;
    int neg_index = mapped_basis[gen_index].second;
    if (pos_index != neg_index) {
      fmpz_set_si(fmpz_mat_entry(boundary_map_matrix, pos_index, gen_index), 1);
      fmpz_set_si(fmpz_mat_entry(boundary_map_matrix, neg_index, gen_index), -1);
    }
  }

  // printf("boundary map matrix:\n");
  // fmpz_mat_print_pretty(boundary_map_matrix);
  // printf("\n");

  // Compute the (right) kernel of the boundary map matrix
  fmpz_mat_t boundary_map_kernel;
  fmpz_mat_init(boundary_map_kernel, b_size, b_size);
  int64_t rank = fmpz_mat_nullspace_mul(boundary_map_kernel, boundary_map_matrix);

  fmpz_mat_clear(boundary_map_matrix);

  fmpz_mat_div_colwise_gcd(boundary_map_kernel);

  fmpz_mat_t trimmed_kernel;
  fmpz_mat_init(trimmed_kernel, b_size, rank);

  for (int row = 0; row < b_size; row++) {
    for (int col = 0; col < rank; col++) {
      fmpz_swap(
        fmpz_mat_entry(trimmed_kernel, row, col),
        fmpz_mat_entry(boundary_map_kernel, row, col)
      );
    }
  }

  fmpz_mat_clear(boundary_map_kernel);

  // printf("boundary map kernel:\n");
  // fmpz_mat_print_pretty(trimmed_kernel);
  // printf("\n");

  DenseBasis out;
  out.set_move(trimmed_kernel);

  return out;
}