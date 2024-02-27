#include "manin_symbol.h"
#include "manin_element.h"
#include "manin_basis.h"

#include <flint/fmpq_mat.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz.h>

#include <iostream>
#include <cassert>

// --- MBEWC functions ---

MBEWC::~MBEWC() {
  fmpq_clear(&coeff);
}

MBEWC MBEWC::negate() const {
  fmpq_t neg_coeff;
  fmpq_init(neg_coeff);
  fmpq_neg(neg_coeff, &(this->coeff));
  return {.basis_index = this->basis_index, .coeff = *neg_coeff};
}

MBEWC MBEWC::scale(const fmpq_t c) const {
  fmpq_t scaled_coeff;
  fmpq_init(scaled_coeff);
  fmpq_mul(scaled_coeff, &(this->coeff), c);
  return {.basis_index = this->basis_index, .coeff = *scaled_coeff};
}

void MBEWC::print() const {
  fmpq_print(&coeff);
  printf(" * [%lld]", basis_index);
}

// --- ManinElement functions ---

ManinElement ManinElement::zero(const int64_t level) {
  std::vector<MBEWC> empty_vec;
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
      fmpq_add(new_coeff, &(it1->coeff), &(it2->coeff));

      if (!fmpq_is_zero(new_coeff)) {
        MBEWC new_MBEWC = {.basis_index = it1->basis_index, .coeff = *new_coeff};
        merged.push_back(new_MBEWC);
      }

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
      fmpq_sub(new_coeff, &(it1->coeff), &(it2->coeff));

      if (!fmpq_is_zero(new_coeff)) {
        MBEWC new_MBEWC = {.basis_index = it1->basis_index, .coeff = *new_coeff};
        merged.push_back(new_MBEWC);
      }

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

  ManinElement result = {.N = N, .components = negated_components};
  result.is_sorted = this->is_sorted;
  return result;
}

ManinElement ManinElement::scale(const fmpq_t c) const {
  std::vector<MBEWC> scaled_components;
  auto scale = [c](MBEWC MBEWC) { return MBEWC.scale(c); };
  std::transform(components.begin(), components.end(), std::back_inserter(scaled_components), scale);

  ManinElement result = {.N = N, .components = scaled_components};
  result.is_sorted = this->is_sorted;
  return result;
}

ManinElement ManinElement::map(std::function<ManinElement(ManinBasisElement)> f, int64_t M) const {
  int64_t out_level = M ? M : N;
  ManinElement result = ManinElement::zero(out_level);

  // printf("ManinElement.map() called with out_level = %lld\n, this = ", out_level);
  // this->print_with_generators();
  // printf("\n");

  std::vector<ManinBasisElement> full_basis = manin_basis(N);
  for (auto MBEWC : components) {
    // [ ] cache result of f?
    auto mbe = full_basis[MBEWC.basis_index];
    auto fmbe = f(mbe);

    // printf("\nin ManinElement::map:\nmbe:");
    // mbe.print_with_indices();
    // printf("\nscale: ");
    // fmpq_print(&(MBEWC.coeff));
    // printf("\nf(mbe):");
    // fmbe.print_with_generators();
    // printf("\n\n");

    result += fmbe.scale(&(MBEWC.coeff));
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
    fmpq_print(&component.coeff);
    printf(" * ");
    manin_basis(N)[component.basis_index].print();
  }
}

std::vector<ManinElement> map_kernel(std::vector<ManinElement> B, std::function<ManinElement(ManinBasisElement)> f, int64_t M) {
  // If the input space is trivial, the result is also trivial.
  if (B.size() == 0) {
    return std::vector<ManinElement>();
  }

  // Otherwise, we find the level from the first element.
  int64_t N = B[0].N;
  std::vector<ManinBasisElement> N_basis = manin_basis(N);
  std::vector<ManinBasisElement> M_basis = manin_basis(M);

  // Construct matrix of B
  // TODO: just make a function that converts between sparse and dense reps of
  // a ManinElement basis
  fmpq_mat_t B_matrix;
  fmpq_mat_init(B_matrix, N_basis.size(), B.size());
  fmpq_mat_zero(B_matrix);

  fmpz_mat_t B_matrix_z;
  fmpz_mat_init(B_matrix_z, N_basis.size(), B.size());
  fmpz_mat_zero(B_matrix_z);

  for (int col = 0; col < B.size(); col++) {
    ManinElement b = B[col];
    // b.print_with_generators();
    // printf("\n");
    for (MBEWC component : b.components) {
      int row = component.basis_index;
      // printf("%d\n", row);
      assert(row < N_basis.size());
      fmpq_set(fmpq_mat_entry(B_matrix, row, col), &component.coeff);
    }
  }

  // XXX: this feels a bit wasteful
  fmpq_mat_get_fmpz_mat_colwise(B_matrix_z, NULL, B_matrix);
  fmpq_mat_clear(B_matrix);

  // printf("basis matrix:\n");
  // fmpz_mat_print_pretty(B_matrix_z);
  // printf("\n");

  // Construct matrix of the map
  fmpq_mat_t map_matrix;
  fmpq_mat_init(map_matrix, M_basis.size(), B.size());
  fmpq_mat_zero(map_matrix);

  fmpz_mat_t map_matrix_z;
  fmpz_mat_init(map_matrix_z, M_basis.size(), B.size());
  fmpz_mat_zero(map_matrix_z);

  for (int col = 0; col < B.size(); col++) {
    // [ ]: maybe inline map() and cache f(mbe)?
    ManinElement fb = B[col].map(f, M);
    // fb.print_with_generators();
    // printf("\n");
    for (MBEWC component : fb.components) {
      int row = component.basis_index;
      // printf("%d\n", row);
      assert(row < M_basis.size());
      fmpq_set(fmpq_mat_entry(map_matrix, row, col), &component.coeff);
    }
  }

  // XXX: this feels a bit wasteful
  fmpq_mat_get_fmpz_mat_rowwise(map_matrix_z, NULL, map_matrix);
  fmpq_mat_clear(map_matrix);

  // printf("map matrix:\n");
  // fmpz_mat_print_pretty(map_matrix_z);
  // printf("\n");

  fmpz_mat_t map_kernel, map_kernel_window;
  fmpz_mat_init(map_kernel, B.size(), B.size());
  int64_t rank = fmpz_mat_nullspace(map_kernel, map_matrix_z);
  fmpz_mat_window_init(map_kernel_window, map_kernel, 0, 0, B.size(), rank);

  fmpz_mat_clear(map_matrix_z);

  // printf("map kernel window:\n");
  // fmpz_mat_print_pretty(map_kernel_window);
  // printf("\n");

  fmpz_mat_t map_kernel_in_orig_basis;
  fmpz_t den;
  fmpz_init(den);
  fmpz_mat_init(map_kernel_in_orig_basis, N_basis.size(), rank);
  fmpz_mat_mul(map_kernel_in_orig_basis, B_matrix_z, map_kernel_window);
  fmpz_mat_content(den, map_kernel_in_orig_basis);
  if (!fmpz_is_zero(den)) {
    fmpz_mat_scalar_divexact_fmpz(map_kernel_in_orig_basis, map_kernel_in_orig_basis, den);
  }

  fmpz_mat_window_clear(map_kernel_window);
  fmpz_mat_clear(map_kernel);

  // printf("map kernel in orig basis:\n");
  // fmpz_mat_print_pretty(map_kernel_in_orig_basis);
  // printf("\n");


  // Convert each column of the kernel to a ManinElement
  std::vector<ManinElement> output;

  // FIXME: multiply result matrix with basis B first.

  for (int col = 0; col < rank; col++) {
    std::vector<MBEWC> components;
    for (int row = 0; row < N_basis.size(); row++) {
      if (!(fmpz_is_zero(fmpz_mat_entry(map_kernel_in_orig_basis, row, col)))) {
        fmpq_t coeff;
        fmpq_init(coeff);
        fmpq_set_fmpz(coeff, fmpz_mat_entry(map_kernel_in_orig_basis, row, col));
        components.push_back({.basis_index = row, .coeff = *coeff});
      }
    }
    ManinElement element = {.N = N, .components = components};
    element.sort();
    output.push_back(element);
  }

  fmpz_mat_clear(map_kernel_in_orig_basis);

  return output;
}