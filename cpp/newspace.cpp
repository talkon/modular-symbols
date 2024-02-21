#include "newspace.h"
#include "modular_symbol.h"
#include "manin_symbol.h"
#include "manin_element.h"

#include <cassert>

ManinElement oldspace_map(ManinGenerator mg, int64_t d, int64_t M) {
  int64_t N = mg.N;
  assert(N % (d * M) == 0);
  IntMatrix2x2 matrix = {.w = d, .x = 0, .y = 0, .z = 1};
  return mg.as_modular_symbol().left_action_by(matrix).to_manin_element(M);
}

