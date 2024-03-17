#include <flint/fmpz_mat.h>
#include <flint/fmpz.h>

#include "debug_utils.h"
#include <stdexcept>

void debug_temp() {

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