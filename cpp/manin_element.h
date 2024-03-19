#ifndef MANIN_ELEMENT_H
#define MANIN_ELEMENT_H

#include <flint/fmpq.h>

#include <functional>
#include <vector>

// Forward declarations
struct ManinBasisElement;

// Manin basis element with coefficient, used as a component in ManinElement.
//
// The generator is stored as the index (in the std::vector returned by ManinGenerator).
// The level N is needed to get an actual ManinGenerator out of an MBEWC.
//
// For now, the coefficient is stored by value to simplify memory management.
struct MBEWC {
  int64_t basis_index;
  fmpq_t coeff;

  MBEWC(int64_t, fmpq_t);

  MBEWC(const MBEWC&);

  ~MBEWC();

  friend auto operator<=> (const MBEWC& left, const MBEWC& right) {
    return left.basis_index <=> right.basis_index;
  }

  // Returns the negation of this MBEWC.
  MBEWC negate() const;

  // Returns this MBEWC, scaled by the given constant.
  MBEWC scale(const fmpq_t) const;

  // Prints this MBEWC.
  void print() const;

  void print_internals() const;
};

// An element of the space of Manin symbols of level N.
// Represented as a sparse linear combination of basis elements.
//
// The indices in `components` must be distinct, and if `is_sorted` is true,
// the indices must be sorted.
//
// XXX: this representation seems inefficient for adding ManinElements;
// perhaps consider other representations? (like a map from index -> coeff)
//
// XXX: actually maybe a dense representation is just faster?
struct ManinElement {
  int64_t N;
  std::vector<MBEWC> components;
  bool is_sorted = false;

  ManinElement(int64_t N, std::vector<MBEWC>);
  ManinElement(const ManinElement&);

  // Zero element of a given level
  static ManinElement zero(const int64_t level);

  // --- Internal representation ---

  // Sorts the components of this element.
  void sort();
  // Marks the components of this element as sorted, without checking.
  void mark_as_sorted_unchecked();

  // --- Overloaded arithmetic operations ---

  ManinElement& operator+= (const ManinElement&);
  ManinElement& operator-= (const ManinElement&);
  friend ManinElement operator+ (const ManinElement&, const ManinElement&);
  friend ManinElement operator- (const ManinElement&, const ManinElement&);

  // --- Mathematical operations ---

  // Returns the negation of this element.
  ManinElement negate() const;

  // Returns this element, scaled by the given constant.
  ManinElement scale(const fmpq_t) const;

  // Returns the result of applying a linear map f to this element.
  // The map f should always return ManinElements of level `M`.
  ManinElement map(std::function<ManinElement(ManinBasisElement)> f, int64_t M) const;

  // --- Displaying ---

  // Prints this Manin element
  void print() const;

  // Prints this Manin element, with index expanded into ManinGenerators
  void print_with_generators() const;

  void print_with_internals() const;
};

#endif // MANIN_ELEMENT_H