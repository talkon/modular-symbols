#include "manin_symbol.h"
#include "manin_element.h"

#include <iostream>
#include <cassert>

MGWC::~MGWC() {
  fmpq_clear(&coeff);
}

MGWC MGWC::negate() const {
  fmpq_t neg_coeff;
  fmpq_init(neg_coeff);
  fmpq_neg(neg_coeff, &(this->coeff));

  MGWC negated = {.index = this->index, .coeff = *neg_coeff};

  return negated;
}

void MGWC::print() const {
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
  assert(is_sorted);
  assert(other.is_sorted);

  std::vector<MGWC> merged;

  auto it1 = this->components.begin();
  auto it2 = other.components.begin();

  while (true) {

    // If we run out of components in `this`, append the remaining part of `other` to `merged`.
    if (it1 == this->components.end()) {
      merged.insert(merged.end(), it2, other.components.end());
      break;
    }

    // If we run out of components in` other`, append the remaining part of `this` to `merged`.
    if (it2 == other.components.end()) {
      merged.insert(merged.end(), it1, this->components.end());
      break;
    }

    if (it1->index < it2->index) {
      merged.push_back(*it1);
      it1++;
    }

    if (it1->index == it2->index) {
      fmpq_t new_coeff;
      fmpq_init(new_coeff);
      fmpq_add(new_coeff, &(it1->coeff), &(it2->coeff));

      if (!fmpq_is_zero(new_coeff)) {
        MGWC new_mgwc = {.index = it1->index, .coeff = *new_coeff};
        merged.push_back(new_mgwc);
      }

      it1++;
      it2++;
    }

    if (it1->index > it2->index) {
      merged.push_back(*it2);
      it2++;
    }
  }

  this->components = merged;
  return *this;
}

ManinElement& ManinElement::operator-= (const ManinElement& other) {
  assert(is_sorted);
  assert(other.is_sorted);

  std::vector<MGWC> merged;

  auto it1 = this->components.begin();
  auto it2 = other.components.begin();

  while (true) {

    // If we run out of components in `this`, negate the remaining part of `other`,
    // and append to `merged`.
    if (it1 == this->components.end()) {
      for (; it2 != other.components.end(); it2++) {
        merged.push_back((*it2).negate());
      }
      break;
    }

    // If we run out of components in` other`, append the remaining part of `this` to `merged`.
    if (it2 == other.components.end()) {
      merged.insert(merged.end(), it1, this->components.end());
      break;
    }

    if (it1->index < it2->index) {
      merged.push_back(*it1);
      it1++;
    } else if (it1->index == it2->index) {
      fmpq_t new_coeff;
      fmpq_init(new_coeff);
      fmpq_sub(new_coeff, &(it1->coeff), &(it2->coeff));

      if (!fmpq_is_zero(new_coeff)) {
        MGWC new_mgwc = {.index = it1->index, .coeff = *new_coeff};
        merged.push_back(new_mgwc);
      }

      it1++;
      it2++;
    } else if (it1->index > it2->index) {
      merged.push_back((*it2).negate());
      it2++;
    }
  }

  this->components = merged;
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