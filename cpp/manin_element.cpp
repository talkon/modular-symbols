#include "manin_symbol.h"
#include "manin_element.h"
#include "manin_basis.h"
#include "debug_utils.h"

#include <flint/fmpq_mat.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz.h>

#include <iostream>
#include <set>
#include <cassert>
#include <stdexcept>


// --- MBEWC functions ---

MBEWC::MBEWC(int64_t bi, fmpq_t c) : basis_index(bi) {
  fmpq_init(coeff);
  fmpq_set(coeff, c);
}

MBEWC::MBEWC(const MBEWC& mbewc) : basis_index(mbewc.basis_index) {
  fmpq_init(coeff);
  fmpq_set(coeff, mbewc.coeff);
}

MBEWC::MBEWC(MBEWC&& mbewc) : basis_index(mbewc.basis_index) {
  fmpq_init(coeff);
  fmpq_swap(coeff, mbewc.coeff);
}

MBEWC& MBEWC::operator=(const MBEWC& mbewc) {
  basis_index = mbewc.basis_index;
  fmpq_init(coeff);
  fmpq_set(coeff, mbewc.coeff);
  return *this;
}

MBEWC& MBEWC::operator=(MBEWC&& mbewc) {
  basis_index = mbewc.basis_index;
  fmpq_init(coeff);
  fmpq_swap(coeff, mbewc.coeff);
  return *this;
}

MBEWC::~MBEWC() {
  fmpq_clear(coeff);
}

MBEWC MBEWC::negate() const {
  fmpq_t neg_coeff;
  fmpq_init(neg_coeff);
  fmpq_neg(neg_coeff, coeff);
  auto out = MBEWC(basis_index, neg_coeff);
  fmpq_clear(neg_coeff);
  return out;
}

MBEWC MBEWC::scale(const fmpq_t c) const {
  fmpq_t scaled_coeff;
  fmpq_init(scaled_coeff);
  fmpq_mul(scaled_coeff, coeff, c);
  auto out = MBEWC(basis_index, scaled_coeff);
  fmpq_clear(scaled_coeff);
  return out;
}

void MBEWC::print() const {
  fmpq_print(coeff);
  printf(" * [%lld]", basis_index);
}

void MBEWC::print_internals() const {
  printf("(");
  uint64_t num = *fmpq_numref(coeff);
  if (COEFF_IS_MPZ(num)) printf("<ref 0x%llx> / ", num);
  else printf("%lld / ", num);

  uint64_t den = *fmpq_denref(coeff);
  if (COEFF_IS_MPZ(den)) printf("<ref 0x%llx>", den);
  else printf("%lld", den);

  printf(") * [%lld]", basis_index);
}

// --- ManinElement functions ---

bool has_duplicate_keys(const std::vector<MBEWC>& components) {
  std::map<uint64_t, int> map;
  for (auto it = components.begin(); it != components.end(); it++) {
    uint64_t key = *fmpq_numref(it->coeff);
    if (COEFF_IS_MPZ(key)) {
      if (auto search = map.find(key); search != map.end()) {
        printf(RED "<warning>" RESET " key %llx duplicated, first index: %d, second index: %lld\n", key, map[key], it->basis_index);
        return true;
      }
      map.insert(std::make_pair(key, it->basis_index));
    }
  }
  return false;
}

ManinElement::ManinElement(int64_t N, std::vector<MBEWC>& components) : N(N) {
  this->components.clear();
  this->components.insert(this->components.end(), components.begin(), components.end());
}

ManinElement::ManinElement(const ManinElement& me)
:
  N(me.N),
  is_sorted(me.is_sorted)
{
  this->components.clear();
  this->components.insert(this->components.end(), me.components.begin(), me.components.end());
}

ManinElement ManinElement::zero(const int64_t level) {
  std::vector<MBEWC> empty_vec;
  ManinElement result = ManinElement(level, empty_vec);
  result.mark_as_sorted_unchecked();
  return result;
}

void ManinElement::mark_as_sorted_unchecked() {
  is_sorted = true;
}

ManinElement& ManinElement::operator+= (const ManinElement& other) {
  assert(is_sorted);
  assert(other.is_sorted);
  assert(N == other.N);

  std::vector<MBEWC> merged;

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

    if (it1->basis_index < it2->basis_index) {
      merged.push_back(*it1);
      it1++;
    } else if (it1->basis_index == it2->basis_index) {
      fmpq_t new_coeff;
      fmpq_init(new_coeff);
      fmpq_add(new_coeff, it1->coeff, it2->coeff);

      if (!fmpq_is_zero(new_coeff)) {
        MBEWC new_MBEWC = MBEWC(it1->basis_index, new_coeff);
        merged.push_back(new_MBEWC);
      }

      fmpq_clear(new_coeff);

      it1++;
      it2++;
    } else if (it1->basis_index > it2->basis_index) {
      auto x = *it2;
      merged.push_back(x);
      it2++;
    }
  }

  // This is fine with move semantics
  this->components = merged;

  return *this;
}

ManinElement& ManinElement::operator-= (const ManinElement& other) {
  assert(is_sorted);
  assert(other.is_sorted);
  assert(N == other.N);

  std::vector<MBEWC> merged;

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

    if (it1->basis_index < it2->basis_index) {
      merged.push_back(*it1);
      it1++;
    } else if (it1->basis_index == it2->basis_index) {
      fmpq_t new_coeff;
      fmpq_init(new_coeff);
      fmpq_sub(new_coeff, it1->coeff, it2->coeff);

      if (!fmpq_is_zero(new_coeff)) {
        MBEWC new_MBEWC = MBEWC(it1->basis_index, new_coeff);
        merged.push_back(new_MBEWC);
      }

      fmpq_clear(new_coeff);

      it1++;
      it2++;
    } else if (it1->basis_index > it2->basis_index) {
      merged.push_back((*it2).negate());
      it2++;
    }
  }

  // This is fine with move semantics
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

ManinElement ManinElement::negate() const {
  std::vector<MBEWC> negated_components;
  auto negate = [](MBEWC MBEWC) { return MBEWC.negate(); };
  std::transform(components.begin(), components.end(), std::back_inserter(negated_components), negate);

  ManinElement result = ManinElement(N, negated_components);
  result.is_sorted = this->is_sorted;
  return result;
}

ManinElement ManinElement::scale(const fmpq_t c) const {
  std::vector<MBEWC> scaled_components;
  auto scale = [c](MBEWC MBEWC) { return MBEWC.scale(c); };
  std::transform(components.begin(), components.end(), std::back_inserter(scaled_components), scale);

  ManinElement result = ManinElement(N, scaled_components);
  result.is_sorted = this->is_sorted;
  return result;
}

ManinElement ManinElement::map(std::function<ManinElement(ManinBasisElement)> f, int64_t M) const {
  int64_t out_level = M ? M : N;
  ManinElement result = ManinElement::zero(out_level);

  auto& full_basis = manin_basis(N);
  for (MBEWC component: components) {
    auto mbe = full_basis[component.basis_index];
    auto fmbe = f(mbe);
    result += fmbe.scale(component.coeff);
  }

  return result;
}

void ManinElement::print() const {
  printf("level: %lld, components:", N);
  for (MBEWC component : components) {
    printf(" + ");
    component.print();
  }
}

void ManinElement::print_with_generators() const {
  for (MBEWC component : components) {
    printf(" + ");
    fmpq_print(component.coeff);
    printf(" * ");
    manin_basis(N)[component.basis_index].print();
  }
}

void ManinElement::print_with_internals() const {
  printf("level: %lld, components:\n", N);
  for (auto it = components.begin(); it != components.end(); it++) {
    printf("  ");
    it->print_internals();
    printf("\n");
  }
}
