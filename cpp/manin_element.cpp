#include "manin_symbol.h"
#include "manin_element.h"
#include "manin_basis.h"

#include <flint/fmpq_mat.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz.h>

#include <iostream>
#include <cassert>

#include "debug_utils.h"

// --- MBEWC functions ---

MBEWC::MBEWC(int64_t bi, fmpq_t c) : basis_index(bi) {
  fmpq_init(coeff);
  fmpq_set(coeff, c);
}

MBEWC::MBEWC(const MBEWC& mbewc) : basis_index(mbewc.basis_index) {
  fmpq_init(coeff);
  fmpq_set(coeff, mbewc.coeff);
}

MBEWC::~MBEWC() {
  // fmpz_t i;
  // fmpz_init_set_ui(i, 2410325428139824902LL);
  // if (fmpq_equal_fmpz(coeff, i)) {
  //   printf("<----------- destroyed!! ----------->\n");
  // }
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

// --- ManinElement functions ---

ManinElement::ManinElement(int64_t N, std::vector<MBEWC> components) : N(N), components(components) {}

ManinElement::ManinElement(const ManinElement& me)
:
  N(me.N),
  is_sorted(me.is_sorted),
  components(me.components)
{
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

      // fmpq_clear(new_coeff);

      it1++;
      it2++;
    } else if (it1->basis_index > it2->basis_index) {
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

  // printf("ManinElement.map() called with out_level = %lld\n, this = ", out_level);
  // this->print_with_generators();
  // printf("\n");

  // if (M == 422) {
  //   DEBUG_INFO_PRINT(3, "1\n");
  // }

  std::vector<ManinBasisElement> full_basis = manin_basis(N);
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

    // if (M == 422) {
    //   DEBUG_INFO(3,
    //     {
    //       component.print();
    //       printf("\n");
    //     }
    //   )

    //   fmbe.scale(component.coeff);

    //   DEBUG_INFO(3,
    //     {
    //       result.print();
    //       printf("\n");
    //     }
    //   )
    // }

    // else {
    result += fmbe.scale(component.coeff);
    // }
  }

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
