#include "manin_symbol.h"
#include "manin_element.h"
#include "manin_basis.h"

#include <flint/fmpq_mat.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz.h>

#include <iostream>
#include <set>
#include <cassert>
#include <stdexcept>

#include "debug_utils.h"
#include "debug_temp.h"

// --- MBEWC functions ---

MBEWC::MBEWC(int64_t bi, fmpq_t c) : basis_index(bi) {
  fmpq_init(coeff);
  fmpq_set(coeff, c);

  uint64_t num = *fmpq_numref(coeff);
  if (COEFF_IS_MPZ(num)) printf("+ <ref 0x%llx>\n", num);

  uint64_t den = *fmpq_denref(coeff);
  if (COEFF_IS_MPZ(den)) printf("+ <ref 0x%llx>\n", den);
}

MBEWC::MBEWC(const MBEWC& mbewc) : basis_index(mbewc.basis_index) {
  fmpq_init(coeff);
  fmpq_set(coeff, mbewc.coeff);

  // printf("MBEWC copy!\n");

  uint64_t old_num = *fmpq_numref(mbewc.coeff);
  uint64_t num = *fmpq_numref(coeff);
  if (COEFF_IS_MPZ(num)) {
    printf("c <ref 0x%llx> -> <ref 0x%llx>\n", old_num, num);
    if (old_num == num) {
      printf(RED "<error>" RESET " attempted move to same address");
      throw std::runtime_error("alkjsdfhah;a");
    }
  //   fmpq_print(mbewc.coeff);
  //   printf(" -> ");
  //   fmpq_print(coeff);
  //   printf("\n");
  }

  uint64_t old_den = *fmpq_denref(mbewc.coeff);
  uint64_t den = *fmpq_denref(coeff);
  if (COEFF_IS_MPZ(den)) printf("c <ref 0x%llx> -> <ref 0x%llx>\n", old_den, den);
}

MBEWC::~MBEWC() {
  uint64_t num = *fmpq_numref(coeff);
  if (COEFF_IS_MPZ(num)) printf("~ <ref 0x%llx>\n", num);

  uint64_t den = *fmpq_denref(coeff);
  if (COEFF_IS_MPZ(den)) printf("~ <ref 0x%llx>\n", den);

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
      if (map.contains(key)) {
        printf(RED "<warning>" RESET " key %llx duplicated, first index: %d, second index: %lld\n", key, map[key], it->basis_index);
        return true;
      }
      map.insert(std::make_pair(key, it->basis_index));
    }
  }
  return false;
}

bool has_duplicate_keys(const std::vector<MBEWC>& c1, const std::vector<MBEWC>& c2) {
  std::map<uint64_t, int> map;
  for (auto it = c1.begin(); it != c1.end(); it++) {
    uint64_t key = *fmpq_numref(it->coeff);
    if (COEFF_IS_MPZ(key)) {
      if (map.contains(key)) {
        printf(RED "<warning>" RESET " key %llx duplicated, first index: %d, second index: %lld\n", key, map[key], it->basis_index);
        return true;
      }
      map.insert(std::make_pair(key, it->basis_index));
    }
  }
  for (auto it = c2.begin(); it != c2.end(); it++) {
    uint64_t key = *fmpq_numref(it->coeff);
    if (COEFF_IS_MPZ(key)) {
      if (map.contains(key)) {
        printf(RED "<warning>" RESET " key %llx duplicated, first index: %d, second index: %lld\n", key, map[key], it->basis_index);
        return true;
      }
      map.insert(std::make_pair(key, it->basis_index));
    }
  }
  return false;
}


ManinElement::ManinElement(int64_t N, std::vector<MBEWC> components) : N(N), components(components) {}

ManinElement::ManinElement(const ManinElement& me)
:
  N(me.N),
  is_sorted(me.is_sorted),
  components(me.components)
{
  // for (auto component : me.components) {
  //   components.push_back(component);
  // }
}

ManinElement ManinElement::zero(const int64_t level) {
  std::vector<MBEWC> empty_vec;
  ManinElement result = ManinElement(level, empty_vec);
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

// BUG: there's something wrong in here about lifetimes?? that causes weird behavior
ManinElement& ManinElement::operator+= (const ManinElement& other) {
  assert(is_sorted);
  assert(other.is_sorted);
  assert(N == other.N);

  std::vector<MBEWC> merged;

  auto it1 = this->components.begin();
  auto it2 = other.components.begin();

  bool found_duplicate_keys = false;
  if (!found_duplicate_keys && has_duplicate_keys(other.components)) {
    printf("other has duplicate keys\n");
    found_duplicate_keys = true;
  }
  if (!found_duplicate_keys && has_duplicate_keys(this->components)) {
    printf("this has duplicate keys\n");
    found_duplicate_keys = true;
  }
  if (!found_duplicate_keys && has_duplicate_keys(this->components, other.components)) {
    printf("operands have duplicate keys\n");
    found_duplicate_keys = true;
  }

  // DEBUG_INFO_PRINT(0, "this: \n");
  // this->print_with_internals();
  // DEBUG_INFO_PRINT(0, "other: \n");
  // other.print_with_internals();

  while (true) {

    // if (N == 422) {
    //   printf(YEL "\nmerged:\n" RESET);
    //   for (auto it = merged.begin(); it != merged.end(); it++) {
    //     printf("  ");
    //     it->print_internals();
    //     printf("\n");
    //   }
    //   printf("\n");
    // }

    // If we run out of components in `this`, append the remaining part of `other` to `merged`.
    if (it1 == this->components.end()) {
      merged.insert(merged.end(), it2, other.components.end());
      if (!found_duplicate_keys && has_duplicate_keys(merged)) {
        printf("found at it2 insert\n");
        found_duplicate_keys = true;
        break;
      }
      break;
    }

    // If we run out of components in` other`, append the remaining part of `this` to `merged`.
    if (it2 == other.components.end()) {
      merged.insert(merged.end(), it1, this->components.end());
      if (!found_duplicate_keys && has_duplicate_keys(merged)) {
        printf("found at it1 insert\n");
        found_duplicate_keys = true;
        break;
      }
      break;
    }

    if (it1->basis_index < it2->basis_index) {
      // printf("it1 at bi %lld\n", it1->basis_index);
      auto x = *it1;
      // printf("dereferenced\n");
      // x.print_internals();
      // printf("attempting push_back\n");
      merged.push_back(x);
      // printf("pushed back\n");
      // merged.back().print_internals();
      // printf("\n");
      if (!found_duplicate_keys && has_duplicate_keys(merged)) {
        printf("found at it1 push_back at bi %lld\n", it1->basis_index);
        found_duplicate_keys = true;
        break;
      }
      it1++;
    } else if (it1->basis_index == it2->basis_index) {
      fmpq_t new_coeff;
      fmpq_init(new_coeff);
      fmpq_add(new_coeff, it1->coeff, it2->coeff);

      if (!fmpq_is_zero(new_coeff)) {
        // printf("add at bi %lld\n", it1->basis_index);
        MBEWC new_MBEWC = MBEWC(it1->basis_index, new_coeff);
        // printf("mbewc created\n");
        // new_MBEWC.print_internals();
        // printf("\n");
        merged.push_back(new_MBEWC);
        // printf("mbewc pushed back\n");
        // merged.back().print_internals();
        // printf("\n");
      }

      fmpq_clear(new_coeff);


      if (!found_duplicate_keys && has_duplicate_keys(merged)) {
        printf("found at add at bi %lld\n", it1->basis_index);
        found_duplicate_keys = true;
        break;
      }

      it1++;
      it2++;
    } else if (it1->basis_index > it2->basis_index) {
      // printf("it2 at bi %lld\n", it2->basis_index);
      auto x = *it2;
      // printf("dereferenced\n");
      // x.print_internals();
      // printf("attempting push_back\n");
      merged.push_back(x);
      // printf("pushed back\n");
      // merged.back().print_internals();
      // printf("\n");
      if (!found_duplicate_keys && has_duplicate_keys(merged)) {
        printf("found at it2 push_back at bi %lld\n", it2->basis_index);
        found_duplicate_keys = true;
        break;
      }
      it2++;
    }
  }

  // printf("merged:\n");
  // for (auto it = merged.begin(); it != merged.end(); it++) {
  //   printf("  ");
  //   it->print_internals();
  //   printf("\n");
  // }

  // Force copying of elements from merged to this.components()
  this->components.clear();
  this->components.insert(this->components.end(), merged.begin(), merged.end());

  // DEBUG_INFO_PRINT(0, "this, after merged: \n");
  // this->print_with_internals();

  // if (found_duplicate_keys) {
  //   DEBUG_INFO_PRINT(0, "this: \n");
  //   this->print_with_internals();
  //   DEBUG_INFO_PRINT(0, "other: \n");
  //   other.print_with_internals();
  //   printf("merged:\n");
  //   for (auto it = merged.begin(); it != merged.end(); it++) {
  //     printf("  ");
  //     it->print_internals();
  //     printf("\n");
  //   }
  //   this->components = merged;
  //   DEBUG_INFO_PRINT(0, "this, after merged: \n");
  //   this->print_with_internals();
  //   throw std::runtime_error("Duplicate keys found");
  // }

  // printf("aaa 1111\n");

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

  std::vector<ManinBasisElement> full_basis = manin_basis(N);
  for (MBEWC component: components) {
    // [ ] cache result of f?
    auto mbe = full_basis[component.basis_index];
    auto fmbe = f(mbe);

    // if (M == 422) printf("cbi: %lld\n", component.basis_index);

    result += fmbe.scale(component.coeff);
    // if (M == 422) printf("aaa 2222\n");
  }

  return result;
}

ManinElement ManinElement::map_debug(std::function<ManinElement(ManinBasisElement)> f, int64_t M) const {
  int64_t out_level = M ? M : N;
  ManinElement result = ManinElement::zero(out_level);

  // printf("ManinElement.map() called with out_level = %lld\n, this = ", out_level);
  // this->print_with_generators();
  // printf("\n");

  // if (M == 422) {
  //   DEBUG_INFO_PRINT(3, "1\n");
  // }

  std::vector<ManinBasisElement> full_basis = manin_basis(N);
  int count = 0;
  for (MBEWC component: components) {
    // [ ] cache result of f?
    auto mbe = full_basis[component.basis_index];
    auto fmbe = f(mbe);

    // printf("\nin ManinElement::map:\nmbe:");
    // mbe.print_with_indices();
    // printf("\nscale: ");
    // fmpq_print(&(MBEWC.coeff));
    // printf("\nf(mbe):");
    // fmbe.print_with_generators();
    // printf("\n\n");

    // if (M == 422) printf("cbi: %lld\n", component.basis_index);

    auto elt = fmbe.scale(component.coeff);

    result += elt;

    if (count == 81) {
      probe_fmpz_freelist(100);
    }

    count++;
    printf("count: %d\n", count);
  }

  printf("count: %d\n", count);

  // printf("result:");
  // result.print_with_generators();
  // printf("\n");

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
