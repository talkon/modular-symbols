#include "manin_symbol.h"
#include "manin_element.h"

#include <iostream>
#include <cassert>

void MGWC::print() {
  fmpq_print(&coeff);
  printf(" * [%lld]", index);
}

ManinElement ManinElement::zero(const int64_t level) {
  std::vector<MGWC> empty_vec;
  ManinElement result = {.N = level, .components = empty_vec};
  result.mark_as_sorted_unchecked();
  return result;
}

void ManinElement::sort() {
  std::sort(components.begin(), components.end());
  is_sorted = true;
}

void ManinElement::mark_as_sorted_unchecked() {
  is_sorted = true;
}

ManinElement& ManinElement::operator+= (const ManinElement& other) {
  // TODO: implement this
  assert(is_sorted);
  assert(other.is_sorted);

  return *this;
}

ManinElement& ManinElement::operator-= (const ManinElement& other) {
  // TODO: implement this
  assert(is_sorted);
  assert(other.is_sorted);

  return *this;
}

ManinElement operator+ (const ManinElement& left, const ManinElement& right) {
  ManinElement result = left;
  result += right;
  return result;
}

ManinElement operator- (const ManinElement& left, const ManinElement& right) {
  ManinElement result = left;
  result -= right;
  return result;
}

void ManinElement::print() const {
  printf("level: %lld, components:", N);
  for (MGWC component : components) {
    printf(" + ");
    component.print();
  }
}

void ManinElement::print_with_generators() const {
  for (MGWC component : components) {
    printf(" + ");
    fmpq_print(&component.coeff);
    printf(" * ");
    manin_generators(N)[component.index].print();
  }
}