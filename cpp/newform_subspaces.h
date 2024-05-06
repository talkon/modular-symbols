#ifndef NEWFORM_SUBSPACES_H
#define NEWFORM_SUBSPACES_H

#include "modular_symbol.h"
#include "manin_element.h"
#include "flint_wrappers.h"

#include <vector>
#include <map>

// Forward declarations
struct ManinBasisElement;

struct Subspace {
  std::vector<ManinElement> basis;
  bool is_newform_subspace;

  int64_t level;
  std::vector<int64_t> atkin_lehner_pos;
  std::vector<int64_t> atkin_lehner_neg;
  int trace_depth = 0;

  // NOTE: Making an assumption here that trace form coeffs fit in 64 bits
  std::map<int64_t, int64_t> trace_form;
  std::map<int64_t, FmpqMatrix> hecke_matrices;
  // std::map<int64_t, FmpzPoly> hecke_minpolys;

  Subspace(
    std::vector<ManinElement> basis,
    bool is_newform_subspace,
    int64_t level,
    std::vector<int64_t> atkin_lehner_pos,
    std::vector<int64_t> atkin_lehner_neg
  ) :
    basis(basis),
    is_newform_subspace(is_newform_subspace),
    level(level),
    atkin_lehner_pos(atkin_lehner_pos),
    atkin_lehner_neg(atkin_lehner_neg)
  {}

  int dimension() const;

  // Computes the next coefficient of the trace form, and return trace_depth.
  int compute_next_trace();

  void compute_trace_until(int depth);

  // Prints information about this subspace
  void print() const;
};

// Computes the action of the Hecke operator T_p on the given Manin basis element.
// p must not divide the level N.
ManinElement& hecke_action(ManinBasisElement, int64_t p);

// Computes the matrix of the action of the Hecke operator T_p.
// p must not divide the level.
FmpqMatrix& hecke_matrix(int64_t level, int64_t p);

// Computes the action of the Atkin-Lehner involution w_q on the given Manin basis element.
// q must be a prime power such that q || N.
ManinElement& atkin_lehner_action(ManinBasisElement, int64_t q);

// Computes the newform subspaces of a given level.
std::vector<Subspace> newform_subspaces(int64_t level, bool dimension_only, int trace_depth, bool prime_opt);

// Computes the dimensions of the newform subspaces of a given level.
std::vector<int> newform_subspace_dimensions(int64_t level);

#endif // NEWFORM_SUBSPACES_H