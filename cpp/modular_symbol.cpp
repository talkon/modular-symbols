#include "modular_symbol.h"
#include "manin_symbol.h"
#include "manin_basis.h"


ModularSymbol ModularSymbol::left_action_by(const IntMatrix2x2 mat) {
  return {
    .a = mat.x * a + mat.y * c,
    .b = mat.x * b + mat.y * d,
    .c = mat.z * a + mat.w * c,
    .d = mat.z * b + mat.w * d
  };
}

ManinElement ModularSymbol::to_manin_element(const int64_t level) {
  return fraction_to_manin_element(b, d, level) - fraction_to_manin_element(a, c, level);
}