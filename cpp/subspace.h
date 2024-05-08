#ifndef SUBSPACE_H
#define SUBSPACE_H

#include "manin_element.h"
#include "flint_wrappers.h"

#include <vector>
#include <map>
#include <optional>

struct Subspace {
  std::vector<ManinElement> basis;
  bool is_newform_subspace;

  int64_t level;
  std::vector<int64_t> atkin_lehner_pos;
  std::vector<int64_t> atkin_lehner_neg;
  int trace_depth = 0;

  std::optional<FmpzPoly> hecke_field_poly;

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

#endif // SUBSPACE_H