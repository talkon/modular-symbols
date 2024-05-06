#include "heilbronn.h"
#include "modular_symbol.h"
#include "manin_element.h"
#include "manin_basis.h"
#include "newspace.h"
#include "linalg.h"
#include "utils.h"
#include "debug_utils.h"
#include "cache_decorator.h"
#include "newform_subspaces.h"
#include "fmpz_mat_helpers.h"
#include "flint_wrappers.h"

#include <flint/ulong_extras.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>

#include <set>
#include <cassert>
#include <cmath>

int Subspace::dimension() const {
  return basis.size();
}

void Subspace::print() const {
  printf("%d:", dimension());

  n_factor_t factors;
  n_factor_init(&factors);
  n_factor(&factors, level, 1);

  for (int i = 0; i < factors.num; i++) {
    int64_t p = factors.p[i];
    if (std::find(atkin_lehner_pos.begin(), atkin_lehner_pos.end(), p) != atkin_lehner_pos.end()) {
      printf("+%lld", p);
    } else {
      printf("-%lld", p);
    }
    if (i < factors.num - 1) printf(",");
  }

  printf(":");
  for (int i = 1; i <= trace_depth; i++) {
    printf("%lld,", trace_form.at(i));
  }
}

ManinElement _impl_hecke_action(ManinBasisElement mbe, int64_t p) {
  int64_t level = mbe.N;
  ManinElement result = ManinElement::zero(level);
  // FIXME: it seems like something is wrong with the computation of Merel matrices
  if (level % p == 0) {
    for (auto mat : heilbronn_merel(p)) {
      ManinGenerator mg = mbe.right_action_by(mat).as_generator();
      // XXX: this += is unnecessarily expensive, so we're using hecke_matrix (below) instead
      result += level_and_index_to_basis(level, mg.index);
    }
  } else {
    for (auto mat : heilbronn_cremona(p)) {
      ManinGenerator mg = mbe.right_action_by(mat).as_generator();
      result += level_and_index_to_basis(level, mg.index);
    }
  }
  return result;
}

ManinElement& hecke_action(ManinBasisElement mbe, int64_t p) {
  static CacheDecorator<ManinElement, ManinBasisElement, int64_t> _cache_hecke_action(_impl_hecke_action);
  return _cache_hecke_action(mbe, p);
}

FmpqMatrix _impl_hecke_matrix(int64_t level, int64_t p) {
  assert(level % p != 0);
  auto basis = manin_basis(level);

  fmpq_mat_t map_of_basis;
  // I think this inits to zero
  fmpq_mat_init(map_of_basis, basis.size(), basis.size());

  for (int col = 0; col < basis.size(); col++) {
    for (auto mat : heilbronn_cremona(p)) {
      auto mg = basis[col].right_action_by(mat).as_generator();
      auto me = level_and_index_to_basis(level, mg.index);

      for (auto component : me.components) {
        fmpq_add(
          fmpq_mat_entry(map_of_basis, component.basis_index, col),
          fmpq_mat_entry(map_of_basis, component.basis_index, col),
          component.coeff
        );
      }
    }
  }

  FmpqMatrix mat;
  mat.set_move(map_of_basis);
  return mat;
}

FmpqMatrix& hecke_matrix(int64_t level, int64_t p) {
  static CacheDecorator<FmpqMatrix, int64_t, int64_t> _cache_hecke_matrix(_impl_hecke_matrix);
  return _cache_hecke_matrix(level, p);
}

ManinElement _impl_atkin_lehner_action(ManinBasisElement mbe, int64_t q) {
  int64_t level = mbe.N;
  utils::XgcdResult xgcd = utils::xgcd(q, level / q);
  IntMatrix2x2 matrix = {.x = q, .y = -xgcd.b, .z = level, .w = xgcd.a * q};
  return mbe.as_modular_symbol().left_action_by(matrix).to_manin_element(level);
}

ManinElement& atkin_lehner_action(ManinBasisElement mbe, int64_t q) {
  static CacheDecorator<ManinElement, ManinBasisElement, int64_t> _cache_atkin_lehner_action(_impl_atkin_lehner_action);
  return _cache_atkin_lehner_action(mbe, q);
}

std::vector<Subspace> newform_subspaces(int64_t level, bool dimension_only, int trace_depth, bool prime_opt) {
  std::vector<ManinElement> basis = newspace_basis(level);

  DEBUG_INFO_PRINT(1, "Starting computation of newform subspaces for level %lld\n", level);

  std::vector<Subspace> remaining;
  remaining.emplace_back(
    basis, false, level, std::vector<int64_t>(), std::vector<int64_t>()
  );

  n_factor_t factors;
  n_factor_init(&factors);
  n_factor(&factors, level, 1);

  for (int i = 0; i < factors.num; i++) {
    int64_t q = n_pow(factors.p[i], factors.exp[i]);
    DEBUG_INFO_PRINT(2, "Decomposing spaces using Atkin-Lehner involution for prime %lld\n", (int64_t) factors.p[i]);
    auto f = [q](ManinBasisElement mbe) { return atkin_lehner_action(mbe, q); };
    std::vector<Subspace> updated;

    for (auto& subspace : remaining) {
      auto splitted = split(subspace.basis, f);

      auto added_pos = subspace.atkin_lehner_pos;
      added_pos.push_back(factors.p[i]);

      auto added_neg = subspace.atkin_lehner_neg;
      added_neg.push_back(factors.p[i]);

      if (splitted.pos_space.size() > 0) {
        updated.emplace_back(splitted.pos_space, false, level, added_pos, subspace.atkin_lehner_neg);
      }

      if (splitted.neg_space.size() > 0) {
        updated.emplace_back(splitted.neg_space, false, level, subspace.atkin_lehner_pos, added_neg);
      }
    }

    remaining = updated;
  }

  // Sturm bound = N * \prod_{p|N} (1 + 1/p) * k / 12, here k = 2
  int sturm_bound = level;
  for (int i = 0; i < factors.num; i++) {
    sturm_bound += sturm_bound / factors.p[i];
  }
  sturm_bound /= 6;

  std::vector<Subspace> done;

  n_primes_t prime_iter;
  n_primes_init(prime_iter);
  int64_t prev_p = 0;

  int iter = 0;

  int manin_basis_size = manin_basis(level).size();
  fmpq_mat_t sum_hecke_mat;
  fmpq_mat_init(sum_hecke_mat, manin_basis_size, manin_basis_size);
  FmpqMatrix sum_hecke;
  sum_hecke.set_move(sum_hecke_mat);

  while (remaining.size() > 0) {
    int64_t p = n_primes_next(prime_iter);

    if (level % p == 0) continue;

    if (p > sturm_bound) {
      // Spaces will not split after Sturm bound
      done.insert(done.end(), remaining.begin(), remaining.end());
      remaining.clear();
      break;
    }

    DEBUG_INFO_PRINT(2, "Decomposing spaces using Hecke operator T_%lld\n", p);

    FmpqMatrix& hecke_mat = hecke_matrix(level, p);
    fmpq_mat_add(sum_hecke.mat, sum_hecke.mat, hecke_mat.mat);

    // if (p < 20) continue;

    std::vector<Subspace> new_remaining;
    std::vector<Subspace> special;

    for (auto& subspace : remaining) {

      if (subspace.dimension() == 1) {
        subspace.is_newform_subspace = true;
        done.emplace_back(subspace);
        continue;
      }

      DecomposeResult dr = decompose(subspace.basis, hecke_mat, dimension_only, prime_opt);
      DEBUG_INFO(2,
        {
          printf("dim %zu -> ", subspace.basis.size());
          for (auto basis : dr.done) {
            printf("%zu,", basis.size());
          }
          if (dr.special.size() > 0) {
            printf("[");
            for (auto basis : dr.special) {
              printf("%zu,", basis.size());
            }
            printf("]");
          }
          if (dr.remaining.size() > 0) {
            printf("(");
            for (auto basis : dr.remaining) {
              printf("%zu,", basis.size());
            }
            printf(")");
          }
          printf("\n");
        }
      )

      for (auto& done_basis : dr.done) {
        done.emplace_back(done_basis, true, level, subspace.atkin_lehner_pos, subspace.atkin_lehner_neg);
      }

      for (auto& special_basis : dr.special) {
        if (iter > 0) special.emplace_back(special_basis, true, level, subspace.atkin_lehner_pos, subspace.atkin_lehner_neg);
        else new_remaining.emplace_back(special_basis, true, level, subspace.atkin_lehner_pos, subspace.atkin_lehner_neg);
      }

      for (auto& remaining_basis: dr.remaining) {
        new_remaining.emplace_back(remaining_basis, false, level, subspace.atkin_lehner_pos, subspace.atkin_lehner_neg);
      }
    }
    remaining = new_remaining;
    DEBUG_INFO_PRINT(2, "Completed p = %lld, done = %zu spaces, special = %zu spaces, remaining = %zu spaces\n", p, done.size(), special.size(), remaining.size());

    if (special.size() > 0 && iter > 0) {
      std::vector<Subspace> special_remaining;
      for (auto& subspace : special ) {
        DecomposeResult dr = decompose(subspace.basis, sum_hecke, dimension_only, prime_opt);
        DEBUG_INFO(2,
          {
            printf("dim %zu -> ", subspace.basis.size());
            for (auto basis : dr.done) {
              printf("%zu,", basis.size());
            }
            if (dr.special.size() > 0) {
              printf("[");
              for (auto basis : dr.special) {
                printf("%zu,", basis.size());
              }
              printf("]");
            }
            if (dr.remaining.size() > 0) {
              printf("(");
              for (auto basis : dr.remaining) {
                printf("%zu,", basis.size());
              }
              printf(")");
            }
            printf("\n");
          }
        )

        for (auto& done_basis : dr.done) {
          done.emplace_back(done_basis, true, level, subspace.atkin_lehner_pos, subspace.atkin_lehner_neg);
        }

        for (auto& special_basis : dr.special) {
          remaining.emplace_back(special_basis, true, level, subspace.atkin_lehner_pos, subspace.atkin_lehner_neg);
        }

        for (auto& remaining_basis: dr.remaining) {
          remaining.emplace_back(remaining_basis, false, level, subspace.atkin_lehner_pos, subspace.atkin_lehner_neg);
        }
      }
      DEBUG_INFO_PRINT(2, "Completed p = %lld (special), done = %zu spaces, remaining = %zu spaces\n", p, done.size(), remaining.size());
    }

    prev_p = p;
    iter++;
  }

  n_primes_clear(prime_iter);

  // This will set the first coefficient of trace form for each subspace to the dimension.
  for (auto& subspace : done) {
    subspace.compute_next_trace();
  }

  auto compare = [&done] (int a, int b) {
    for (int i = 1; i <= done.at(a).trace_depth && i <= done.at(b).trace_depth; i++) {
      if (done.at(a).trace_form.at(i) == done.at(b).trace_form.at(i)) continue;
      else return done.at(a).trace_form.at(i) < done.at(b).trace_form.at(i);
    }
    return false;
  };

  int num_subspaces = done.size();
  std::vector<int> subspace_order(num_subspaces);
  for (int i = 0; i < num_subspaces; i++) {
    subspace_order[i] = i;
  }

  std::sort(subspace_order.begin(), subspace_order.end(), compare);

  DEBUG_INFO_PRINT(1, "Starting computation of trace forms for level %lld\n", level);

  if (!dimension_only) {
    int next_depth = 1;
    while (true) {
      next_depth++;
      std::set<int> next_trace_needed;

      for (int i = 0; i < num_subspaces - 1; i++) {
        if (done.at(subspace_order[i]).trace_form == done.at(subspace_order[i+1]).trace_form) {
          next_trace_needed.insert(i);
          next_trace_needed.insert(i+1);
        }
      }

      if (next_trace_needed.empty()) break;

      DEBUG_INFO_PRINT(2,
        "Trace depth: %d, subspaces remaining: %zu\n",
        next_depth,
        next_trace_needed.size()
      )

      for (int i : next_trace_needed) {
        DEBUG_INFO_PRINT(5,
          " Subspace %d, dimension %d, depth %d\n",
          subspace_order[i],
          done.at(subspace_order[i]).dimension(),
          done.at(subspace_order[i]).trace_depth + 1
        );
        done.at(subspace_order[i]).compute_next_trace();
      }

      // This is a bit inefficient, but this is not taking up much execution time anyways.
      std::sort(subspace_order.begin(), subspace_order.end(), compare);
    }

    DEBUG_INFO_PRINT(2, "Computing remaining traces up to depth %d\n", trace_depth);

    for (auto& subspace : done) {
        DEBUG_INFO_PRINT(3,
          " Subspace dimension %d, until depth %d\n",
          subspace.dimension(),
          trace_depth
        );
      subspace.compute_trace_until(trace_depth);
    }
  }

  std::vector<Subspace> sorted;
  for (auto i : subspace_order) {
    sorted.emplace_back(done.at(i));
  }

  return sorted;
}

std::vector<int> newform_subspace_dimensions(int64_t level) {
  auto nss = newform_subspaces(level, true, 0, true);
  std::vector<int> sizes;
  for (auto ns : nss) {
    sizes.push_back(ns.dimension());
  }
  std::sort(sizes.begin(), sizes.end());
  return sizes;
}

int Subspace::compute_next_trace() {
  int n = trace_depth + 1;
  int dim = dimension();
  assert(dim > 0);
  int64_t level = basis[0].N;

  n_factor_t factors;
  n_factor_init(&factors);
  n_factor(&factors, n, 1);

  int g = n_gcd(n, level);

  if (g == 1) {
    fmpq_mat_t f_matrix;
    fmpq_mat_init(f_matrix, dim, dim);

    if (n == 1) {
      fmpq_mat_one(f_matrix);
    }
    else if (n_is_prime(n)) {
      int p = n;
      // TODO: There's a lot of copied code here -- might be possible to refactor
      std::vector<ManinBasisElement> N_basis = manin_basis(level);

      // Construct matrix of B
      fmpq_mat_t B_matrix;
      fmpq_mat_init(B_matrix, N_basis.size(), basis.size());
      fmpq_mat_zero(B_matrix);

      fmpz_mat_t B_matrix_z;
      fmpz_mat_init(B_matrix_z, N_basis.size(), basis.size());
      fmpz_mat_zero(B_matrix_z);

      for (int col = 0; col < basis.size(); col++) {
        ManinElement b = basis[col];
        for (MBEWC component : b.components) {
          int row = component.basis_index;
          assert(row < N_basis.size());
          fmpq_set(fmpq_mat_entry(B_matrix, row, col), component.coeff);
        }
      }

      fmpq_mat_get_fmpz_mat_colwise(B_matrix_z, NULL, B_matrix);
      fmpq_mat_clear(B_matrix);

      DEBUG_INFO(6,
        {
          printf("B_matrix_z: ");
          fmpz_mat_print_dimensions(B_matrix_z);
          printf("\n");
        }
      )

      // Finds pivot rows of the matrix B.
      std::vector<int> pivots;
      fmpz* pivot_coeffs = _fmpz_vec_init(basis.size());

      int current_col = 0;
      for (int row = 0; row < N_basis.size(); row++) {
        bool is_pivot = true;
        for (int col = 0; col < basis.size(); col++) {
          if (col != current_col && !fmpz_is_zero(fmpz_mat_entry(B_matrix_z, row, col))) {
            is_pivot = false;
            break;
          }
        }

        if (!is_pivot) continue;

        if (fmpz_is_zero(fmpz_mat_entry(B_matrix_z, row, current_col))) continue;

        pivots.push_back(row);
        fmpz_set(pivot_coeffs + current_col, fmpz_mat_entry(B_matrix_z, row, current_col));
        current_col++;

        if (current_col == basis.size()) break;
      }

      // Construct matrix of the linear map f acting on B

      FmpqMatrix& hecke_mat = hecke_matrix(level, p);

      fmpq_mat_t pivot_rows;
      fmpq_mat_init(pivot_rows, basis.size(), N_basis.size());

      for (int row = 0; row < basis.size(); row++) {
        for (int col = 0; col < N_basis.size(); col++) {
          fmpq_div_fmpz(
            fmpq_mat_entry(pivot_rows, row, col),
            fmpq_mat_entry(hecke_mat.mat, pivots[row], col),
            (pivot_coeffs + row)
          );
        }
      }

      DEBUG_INFO(6,
        {
          printf("pivot_rows: ");
          fmpq_mat_print_dimensions(pivot_rows);
          printf("\n");
        }
      )

      fmpq_mat_mul_fmpz_mat(f_matrix, pivot_rows, B_matrix_z);
      fmpq_mat_clear(pivot_rows);

      DEBUG_INFO(6,
        {
          printf("f_matrix: ");
          fmpq_mat_print_dimensions(f_matrix);
          printf("\n");
        }
      )

      _fmpz_vec_clear(pivot_coeffs, basis.size());
    }
    else if (factors.num == 1) {
      int64_t p = factors.p[0];
      int exp = factors.exp[0];
      int n_over_p = n_pow(p, exp - 1);
      fmpq_mat_mul(f_matrix, hecke_matrices.at(n_over_p).mat, hecke_matrices.at(p).mat);
      if (level % p != 0) {
        fmpq_mat_t temp;
        fmpz_t P;
        fmpq_mat_init(temp, dim, dim);
        fmpz_init_set_si(P, p);

        fmpq_mat_scalar_mul_fmpz(temp, hecke_matrices.at(n_over_p / p).mat, P);
        fmpq_mat_sub(f_matrix, f_matrix, temp);

        fmpq_mat_clear(temp);
        fmpz_clear(P);
      }
    }
    else {
      int64_t q = n_pow(factors.p[0], factors.exp[0]);
      fmpq_mat_mul(f_matrix, hecke_matrices.at(n / q).mat, hecke_matrices.at(q).mat);
    }
    fmpq_t trace;
    fmpq_init(trace);
    fmpq_mat_trace(trace, f_matrix);
    assert(fmpz_is_one(fmpq_denref(trace)));

    int trace_int = fmpz_get_si(fmpq_numref(trace));
    trace_form.insert(std::make_pair(n, trace_int));
    fmpq_clear(trace);

    FmpqMatrix hecke_matrix;
    hecke_matrix.set_move(f_matrix);
    hecke_matrices.insert(std::make_pair(n, hecke_matrix));

    trace_depth++;
  }
  else {
    // For primes p | N:
    // - if p^2 | N, then the action of the Hecke operator T_p is zero, and
    // - if p^2 \nmid N, then the action of the Hecke operator T_p is -W_p, where W_p is the Atkin-Lehner involution.
    // See Prop 13.3.4 in (Cohen, Stromberg)
    int64_t x = 1;
    int sign = 1;

    for (int i = 0; i < factors.num; i++) {
      int64_t p = factors.p[i];
      int exp = factors.exp[i];

      if (level % p == 0) {
        if (level % (p * p) == 0) sign = 0;
        else if (
          exp % 2 == 1
          && std::find(atkin_lehner_pos.begin(), atkin_lehner_pos.end(), p) != atkin_lehner_pos.end()
        ) {
          sign = -sign;
        }
      }
      else {
        x *= n_pow(p, exp);
      }
    }

    int trace_int = trace_form.at(x) * sign;
    trace_form.insert(std::make_pair(n, trace_int));
    trace_depth++;
  }

  return trace_depth;
}

void Subspace::compute_trace_until(int depth) {
  while (trace_depth < depth) compute_next_trace();
}
