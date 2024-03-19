#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>

#include <stdexcept>
#include <cassert>
#include <map>

#include "debug_utils.h"
#include "manin_element.h"

void debug_temp() {
  fmpq_t large, to_be_freed;
  fmpq_init(large);
  fmpq_init(to_be_freed);
  fmpq_set_str(large, "100000000000000000000", 10);
  std::vector<MBEWC> vec;

  fmpq_set(to_be_freed, large);
  for (int i = 0; i < 10000; i++)
    fmpq_clear(to_be_freed);

  for (int i = 0; i < 33; i++) {
    fmpq_t x;
    fmpq_init(x);
    fmpq_add_si(x, large, i);
    auto y = MBEWC(0, x);
    vec.push_back(y);
    printf(YEL "\nvec:\n" RESET);
    for (auto it = vec.begin(); it != vec.end(); it++) {
      printf("  ");
      it->print_internals();
      printf("\n");
    }
    printf("\n");
  }

  fmpq_clear(large);
}

void probe_fmpz_freelist(int depth) {
  fmpz_t large;
  fmpz_set_str(large, "100000000000000000000", 10);

  fmpz* vec = _fmpz_vec_init(depth);
  std::map<uint64_t, int> map;

  for (int i = 0; i < depth; i++) {
    fmpz_set(vec + i, large);
    uint64_t val = *(vec + i);
    assert(COEFF_IS_MPZ(val));
    if (map.contains(val)) {
      printf(RED "<error>" RESET " key %llx duplicated, first index: %d, second index: %d\n", val, map[val], i);
      throw std::runtime_error("duplicates found");
    }
    map.insert(std::make_pair(val, i));
  }

  _fmpz_vec_clear(vec, depth);

  printf(GRN "<info>" RESET " no duplicates found to depth %d\n", depth);
}

void check_status() {

  // printf(YEL "<starting debug_temp>\n" RESET);

  fmpz_mat_t m, k, window;
  fmpz_mat_init(m, 3, 4);
  fmpz_mat_init(k, 4, 4);
  fmpz_mat_zero(m);

  fmpz_set_str(fmpz_mat_entry(m, 0, 0), "981911205996571756", 10);
  fmpz_set_str(fmpz_mat_entry(m, 0, 1), "240491096708991556", 10);
  fmpz_set_str(fmpz_mat_entry(m, 0, 2), "-16813084898157095", 10);
  fmpz_set_str(fmpz_mat_entry(m, 0, 3), "1189660152027260618", 10);
  fmpz_set_str(fmpz_mat_entry(m, 1, 0), "-566527823500001314", 10);
  fmpz_set_str(fmpz_mat_entry(m, 1, 1), "-131592927616542154", 10);
  fmpz_set_str(fmpz_mat_entry(m, 1, 2), "9889030447594245", 10);
  fmpz_set_str(fmpz_mat_entry(m, 1, 3), "-686023046642061734", 10);
  fmpz_set_str(fmpz_mat_entry(m, 2, 0), "1305348122731461211", 10);
  fmpz_set_str(fmpz_mat_entry(m, 2, 1), "311297869126224287", 10);
  fmpz_set_str(fmpz_mat_entry(m, 2, 2), "-22198802700696000", 10);
  fmpz_set_str(fmpz_mat_entry(m, 2, 3), "1608323242553626377", 10);

  // fmpz_mat_print_pretty(m);
  int rank = fmpz_mat_nullspace(k, m);
  // printf("\nrank: %d\n", rank);

  fmpz_t den;

  for (int col = 0; col < rank; col++) {
    fmpz_mat_window_init(window, k, 0, col, 4, col+1);
    fmpz_mat_content(den, window);
    if (!fmpz_is_zero(den)) {
      fmpz_mat_scalar_divexact_fmpz(window, window, den);
    }
    fmpz_mat_window_clear(window);
  }

  // fmpz_mat_print_pretty(k);
  // printf("\n");

  if (fmpz_equal_si(fmpz_mat_entry(k, 0, 0), 90351233434514035)) {
    printf(GRN "<correct>\n" RESET);
  } else {
    printf(RED "<incorrect>\n" RESET);
  }

  fmpz_mat_clear(m);
  fmpz_mat_clear(k);
  fmpz_clear(den);
}