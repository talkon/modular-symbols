#ifndef NEWFORM_SUBSPACES_H
#define NEWFORM_SUBSPACES_H

#include "modular_symbol.h"
#include "manin_element.h"

#include <vector>

// Forward declarations
struct ManinBasisElement;

struct Subspace {
  std::vector<ManinElement> basis;
  bool is_newform_subspace;

  int64_t level;
  std::vector<int64_t> atkin_lehner_pos;
  std::vector<int64_t> atkin_lehner_neg;
  std::vector<int64_t> trace_form;
  // TODO: Hecke action polynomials + matrices? need to figure out how to store this

  int dimension();

  Subspace(
    std::vector<ManinElement> basis,
    bool is_newform_subspace,
    int64_t level,
    std::vector<int64_t> atkin_lehner_pos,
    std::vector<int64_t> atkin_lehner_neg,
    std::vector<int64_t> trace_form
  ) :
    basis(basis),
    is_newform_subspace(is_newform_subspace),
    level(level),
    atkin_lehner_pos(atkin_lehner_pos),
    atkin_lehner_neg(atkin_lehner_neg),
    trace_form(trace_form)
  {}
};

// Computes the action of the Hecke operator T_p on the given Manin basis element.
// p must not divide the level N.
ManinElement hecke_action(ManinBasisElement, int64_t p);

// Computes the action of the Atkin-Lehner involution w_q on the given Manin basis element.
// q must be a prime power such that q || N.
ManinElement atkin_lehner_action(ManinBasisElement, int64_t q);

// Computes the newform subspaces of a given level.
std::vector<Subspace> newform_subspaces(int64_t level);

// Computes the dimensions of the newform subspaces of a given level.
std::vector<int> newform_subspace_dimensions(int64_t level);

#endif // NEWFORM_SUBSPACES_H