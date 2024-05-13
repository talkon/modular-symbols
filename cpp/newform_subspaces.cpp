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
#include "subspace.h"

#include <flint/ulong_extras.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>

#include <set>
#include <cassert>
#include <cmath>

ManinElement _impl_hecke_action(ManinBasisElement mbe, int64_t p) {
  int64_t level = mbe.N;
  ManinElement result = ManinElement::zero(level);
  // FIXME: it seems like something is wrong with the computation of Merel matrices
  if (level % p == 0) {
    for (auto mat : heilbronn_merel(p)) {
      ManinGenerator mg = mbe.right_action_by(mat).as_generator();
      result += level_and_index_to_basis(level, mg.index);
    }
  } else {
    for (auto mat : heilbronn_cremona(p)) {
      ManinGenerator mg = mbe.right_action_by(mat).as_generator();
      // XXX: this += is unnecessarily expensive, so we're using hecke_matrix (below) instead
      result += level_and_index_to_basis(level, mg.index);
    }
  }
  return result;
}

ManinElement& hecke_action(ManinBasisElement mbe, int64_t p) {
  static CacheDecorator<ManinElement, ManinBasisElement, int64_t> _cache_hecke_action(_impl_hecke_action);
  return _cache_hecke_action(mbe, p);
}

FmpqMatrix hecke_matrix(int64_t level, int64_t p) {
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

// FmpqMatrix& hecke_matrix(int64_t level, int64_t p) {
//   static CacheDecorator<FmpqMatrix, int64_t, int64_t> _cache_hecke_matrix(_impl_hecke_matrix);
//   return _cache_hecke_matrix(level, p);
// }

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

std::vector<Subspace> newform_subspaces(int64_t level, bool dimension_only, int min_trace_depth, int max_trace_depth, bool prime_opt) {

  if (level <= 10) return std::vector<Subspace>();

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

    FmpqMatrix hecke_mat = hecke_matrix(level, p);
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

      DecomposeResult dr = decompose(subspace, hecke_mat, dimension_only, prime_opt);
      DEBUG_INFO(2,
        {
          printf("dim %zu -> ", subspace.basis.size());
          for (auto& space : dr.done) {
            printf("%zu,", space.basis.size());
          }
          if (dr.special.size() > 0) {
            printf("[");
            for (auto& space : dr.special) {
              printf("%zu,", space.basis.size());
            }
            printf("]");
          }
          if (dr.remaining.size() > 0) {
            printf("(");
            for (auto& space : dr.remaining) {
              printf("%zu,", space.basis.size());
            }
            printf(")");
          }
          printf("\n");
        }
      )

      for (auto& space : dr.done) {
        done.emplace_back(space);
      }

      for (auto& space : dr.special) {
        if (iter > 0) special.emplace_back(space);
        else new_remaining.emplace_back(space);
      }

      for (auto& space: dr.remaining) {
        new_remaining.emplace_back(space);
      }
    }
    remaining = new_remaining;
    DEBUG_INFO_PRINT(2, "Completed p = %lld, done = %zu spaces, special = %zu spaces, remaining = %zu spaces\n", p, done.size(), special.size(), remaining.size());

    if (special.size() > 0 && iter > 0) {
      std::vector<Subspace> special_remaining;
      for (auto& subspace : special ) {
        DecomposeResult dr = decompose(subspace, sum_hecke, dimension_only, prime_opt);
        DEBUG_INFO(2,
          {
            printf("dim %zu -> ", subspace.basis.size());
            for (auto& space : dr.done) {
              printf("%zu,", space.basis.size());
            }
            if (dr.special.size() > 0) {
              printf("[");
              for (auto& space : dr.special) {
                printf("%zu,", space.basis.size());
              }
              printf("]");
            }
            if (dr.remaining.size() > 0) {
              printf("(");
              for (auto& space : dr.remaining) {
                printf("%zu,", space.basis.size());
              }
              printf(")");
            }
            printf("\n");
          }
        )

        for (auto& space : dr.done) {
          done.emplace_back(space);
        }

        for (auto& space : dr.special) {
          remaining.emplace_back(space);
        }

        for (auto& space: dr.remaining) {
          remaining.emplace_back(space);
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
    subspace.set_first_trace();
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
    int next_depth = 2;

    while (max_trace_depth == -1 || next_depth <= max_trace_depth) {

      std::vector<bool> next_trace_needed(num_subspaces);
      for (int i = 0; i < num_subspaces; i++) {
        next_trace_needed[i] = false;
      }

      bool hecke_mat_needed = next_depth <= min_trace_depth;

      for (int i = 0; i < num_subspaces - 1; i++) {
        if (done.at(subspace_order[i]).trace_form == done.at(subspace_order[i+1]).trace_form) {
          next_trace_needed[i] = true;
          next_trace_needed[i+1] = true;
          hecke_mat_needed = true;
        }
      }

      if (!hecke_mat_needed) break;

      int count = 0;
      for (int i = 0; i < num_subspaces; i++) {
        if(next_trace_needed[i]) count++;
      }

      if (next_depth <= min_trace_depth) {
        DEBUG_INFO_PRINT(2,
          "Trace depth: %d, for all subspaces (%zu subspaces)\n",
          next_depth,
          num_subspaces
        )
      } else {
        DEBUG_INFO_PRINT(2,
          "Trace depth: %d, subspaces remaining: %zu\n",
          next_depth,
          count
        )
      }

      FmpqMatrix hecke_mat;
      if (n_is_prime(next_depth) && n_gcd(level, next_depth) == 1) {
        hecke_mat = hecke_matrix(level, next_depth);
      }

      for (int i = 0; i < num_subspaces; i++) {
        if (next_trace_needed[i]) {
          DEBUG_INFO_PRINT(5,
            " Subspace %d, dimension %d, depth %d\n",
            subspace_order[i],
            done.at(subspace_order[i]).dimension(),
            next_depth
          );
          done.at(subspace_order[i]).next_trace(next_depth, hecke_mat, max_trace_depth);
        } else if (next_depth <= min_trace_depth) {
          DEBUG_INFO_PRINT(5,
            " Subspace %d, dimension %d, depth %d\n",
            subspace_order[i],
            done.at(subspace_order[i]).dimension(),
            next_depth
          );
          done.at(subspace_order[i]).next_trace(next_depth, hecke_mat, min_trace_depth);
          if (next_depth == min_trace_depth) {
            done.at(subspace_order[i]).clear_hecke_matrices();
          }
        } else {
          done.at(subspace_order[i]).clear_hecke_matrices();
        }
      }
      // This is a bit inefficient, but this is not taking up much execution time anyways.
      std::sort(subspace_order.begin(), subspace_order.end(), compare);

      next_depth++;
    }
  }

  std::vector<Subspace> sorted;
  for (auto i : subspace_order) {
    sorted.emplace_back(done.at(i));
  }

  return sorted;
}

std::vector<int> newform_subspace_dimensions(int64_t level) {
  auto nss = newform_subspaces(level, true, 0, 0, true);
  std::vector<int> sizes;
  for (auto ns : nss) {
    sizes.push_back(ns.dimension());
  }
  std::sort(sizes.begin(), sizes.end());
  return sizes;
}
