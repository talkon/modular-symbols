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

#include <flint/ulong_extras.h>
#include <cmath>

int Subspace::dimension() {
  return basis.size();
}

// TODO: consider caching this, or really, just caching its matrix.

ManinElement _impl_hecke_action(ManinBasisElement mbe, int64_t p) {
  int64_t level = mbe.N;
  ManinElement result = ManinElement::zero(level);
  for (auto mat : heilbronn_cremona(p)) {
    ManinGenerator mg = mbe.right_action_by(mat).as_generator();
    result += level_and_index_to_basis(level, mg.index);
  }
  return result;
}

ManinElement hecke_action(ManinBasisElement mbe, int64_t p) {
  static CacheDecorator<ManinElement, ManinBasisElement, int64_t> _cache_hecke_action(_impl_hecke_action);
  return _cache_hecke_action(mbe, p);
}

ManinElement _impl_atkin_lehner_action(ManinBasisElement mbe, int64_t q) {
  int64_t level = mbe.N;
  utils::XgcdResult xgcd = utils::xgcd(q, level / q);
  IntMatrix2x2 matrix = {.x = q, .y = -xgcd.b, .z = level, .w = xgcd.a * q};
  return mbe.as_modular_symbol().left_action_by(matrix).to_manin_element(level);
}

ManinElement atkin_lehner_action(ManinBasisElement mbe, int64_t q) {
  static CacheDecorator<ManinElement, ManinBasisElement, int64_t> _cache_atkin_lehner_action(_impl_atkin_lehner_action);
  return _cache_atkin_lehner_action(mbe, q);
}

std::vector<Subspace> newform_subspaces(int64_t level) {
  std::vector<ManinElement> basis = newspace_basis(level);

  DEBUG_INFO_PRINT(1, "Starting computation of newform subspaces for level %lld\n", level);

  std::vector<Subspace> remaining;
  remaining.emplace_back(
    basis, false, level, std::vector<int64_t>(), std::vector<int64_t>(), std::vector<int64_t>()
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
        updated.emplace_back(splitted.pos_space, false, level, added_pos, subspace.atkin_lehner_neg, subspace.trace_form);
      }

      if (splitted.neg_space.size() > 0) {
        updated.emplace_back(splitted.neg_space, false, level, subspace.atkin_lehner_pos, added_neg, subspace.trace_form);
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

  while (remaining.size() > 0) {
    int64_t p = n_primes_next(prime_iter);

    // TODO:
    if (level % p == 0) continue;

    if (p > sturm_bound) {
      // Spaces will not split after Sturm bound
      done.insert(done.end(), remaining.begin(), remaining.end());
      remaining.clear();
      break;
    }

    if (prev_p) {
      DEBUG_INFO_PRINT(2, "Decomposing spaces using Hecke operators T_%lld + T_%lld\n", prev_p, p);
    } else {
      DEBUG_INFO_PRINT(2, "Decomposing spaces using Hecke operators T_%lld\n", p);
    }

    auto f = [prev_p, p](ManinBasisElement mbe) {
      if (prev_p) {
        return hecke_action(mbe, p) + hecke_action(mbe, prev_p);
      } else {
        return hecke_action(mbe, p);
      }
    };
    std::vector<Subspace> new_remaining;
    // XXX: This causes the action of `f` to be recomputed many times.
    for (auto& subspace : remaining) {

      if (subspace.dimension() == 1) {
        subspace.is_newform_subspace = true;
        done.emplace_back(subspace);
        continue;
      }

      DecomposeResult dr = decompose(subspace.basis, f);
      DEBUG_INFO(2,
        {
          printf("dim %zu -> ", subspace.basis.size());
          for (auto basis : dr.done) {
            printf("%zu,", basis.size());
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
        done.emplace_back(done_basis, true, level, subspace.atkin_lehner_pos, subspace.atkin_lehner_neg, std::vector<int64_t>());
      }

      for (auto& remaining_basis: dr.remaining) {
        new_remaining.emplace_back(remaining_basis, false, level, subspace.atkin_lehner_pos, subspace.atkin_lehner_neg, std::vector<int64_t>());
      }
    }
    remaining = new_remaining;
    DEBUG_INFO_PRINT(2, "Completed p = %lld, done = %zu spaces, remaining = %zu spaces\n", p, done.size(), remaining.size());

    prev_p = p;
  }

  n_primes_clear(prime_iter);

  return done;
}

std::vector<int> newform_subspace_dimensions(int64_t level) {
  auto nss = newform_subspaces(level);
  std::vector<int> sizes;
  for (auto ns : nss) {
    sizes.push_back(ns.dimension());
  }
  std::sort(sizes.begin(), sizes.end());
  return sizes;
}
