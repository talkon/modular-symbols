#include "linalg.h"
#include "manin_basis.h"
#include "manin_element.h"
#include "debug_utils.h"
#include "newform_subspaces.h"
#include "fmpz_mat_helpers.h"
#include "subspace_basis.h"

#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/ulong_extras.h>

#include <vector>
#include <cassert>
#include <stdexcept>

DenseBasis map_kernel(DenseBasis& B, std::function<ManinElement(ManinBasisElement)> f, int64_t N, int64_t M) {
  // If the input space is trivial, the result is also trivial.
  int b_size = B.mat->c;
  if (b_size == 0) {
    return B;
  }

  std::vector<ManinBasisElement> N_basis = manin_basis(N);
  std::vector<ManinBasisElement> M_basis = manin_basis(M);

  // Construct matrix of the map
  fmpq_mat_t map_matrix;
  fmpq_mat_init(map_matrix, M_basis.size(), b_size);
  fmpq_mat_zero(map_matrix);

  fmpq_mat_t map_of_basis;
  fmpq_mat_init(map_of_basis, M_basis.size(), N_basis.size());
  for (int col = 0; col < N_basis.size(); col++) {
    ManinElement fb = f(N_basis[col]);

    for (MBEWC component : fb.components) {
      int row = component.basis_index;
      assert(row < M_basis.size());
      fmpq_set(fmpq_mat_entry(map_of_basis, row, col), component.coeff);
    }
  }

  DEBUG_INFO(4,
    {
      printf("map_of_basis: ");
      fmpq_mat_print_dimensions(map_of_basis);
      printf("\n");
    }
  )

  fmpq_mat_mul_fmpz_mat(map_matrix, map_of_basis, B.mat);
  fmpq_mat_clear(map_of_basis);

  fmpz_mat_t map_matrix_z;
  fmpz_mat_init(map_matrix_z, M_basis.size(), b_size);
  fmpz_mat_zero(map_matrix_z);

  // XXX: the conversions between fmpz and fmpq matrices might use too much memory?
  fmpq_mat_get_fmpz_mat_rowwise(map_matrix_z, NULL, map_matrix);
  fmpq_mat_clear(map_matrix);

  fmpz_mat_t map_kernel, map_kernel_window;
  fmpz_mat_init(map_kernel, b_size, b_size);

  DEBUG_INFO(4,
    {
      printf("map_matrix_z: ");
      fmpz_mat_print_dimensions(map_matrix_z);
      printf("\n");
    }
  )

  int64_t rank = fmpz_mat_nullspace_mul(map_kernel, map_matrix_z);

  fmpz_mat_clear(map_matrix_z);

  fmpz_mat_window_init(map_kernel_window, map_kernel, 0, 0, b_size, rank);

  DEBUG_INFO(4,
    {
      printf("kernel: ");
      fmpz_mat_print_dimensions(map_kernel_window);
      printf("\n");
    }
  )

  fmpz_mat_div_colwise_gcd(map_kernel_window);

  fmpz_mat_t map_kernel_in_orig_basis;
  fmpz_mat_init(map_kernel_in_orig_basis, N_basis.size(), rank);

  DEBUG_INFO(4,
    {
      printf("kernel (cleared): ");
      fmpz_mat_print_dimensions(map_kernel_window);
      printf("\n");
    }
  )

  // XXX: This multiplication is still slow and uses too much memory.
  fmpz_mat_mul(map_kernel_in_orig_basis, B.mat, map_kernel_window);

  DEBUG_INFO(4,
    {
      printf("result: ");
      fmpz_mat_print_dimensions(map_kernel_in_orig_basis);
      printf("\n");
    }
  )

  B.clear();

  fmpz_mat_window_clear(map_kernel_window);
  fmpz_mat_clear(map_kernel);

  fmpz_mat_div_colwise_gcd(map_kernel_in_orig_basis);

  DEBUG_INFO(4,
    {
      printf("result (cleared): ");
      fmpz_mat_print_dimensions(map_kernel_in_orig_basis);
      printf("\n");
    }
  )

  FmpzMatrix out;
  out.set_move(map_kernel_in_orig_basis);

  return out;
}

SplitResult SplitResult::empty(int64_t level) {
  return {
    .pos_space = dense_empty(level),
    .neg_space = dense_empty(level)
  };
}

SplitResult split(DenseBasis& B, std::function<ManinElement (ManinBasisElement)> f, int64_t N) {
  // If the input space is trivial, the result is also trivial.
  int b_size = B.mat->c;
  if (b_size == 0) {
    return SplitResult::empty(N);
  }

  std::vector<ManinBasisElement> N_basis = manin_basis(N);

  // Finds pivot rows of the matrix B.
  std::vector<int> pivots;
  fmpz* pivot_coeffs = _fmpz_vec_init(b_size);

  int current_col = 0;
  for (int row = 0; row < N_basis.size(); row++) {
    bool is_pivot = true;
    for (int col = 0; col < b_size; col++) {
      if (col != current_col && !fmpz_is_zero(fmpz_mat_entry(B.mat, row, col))) {
        is_pivot = false;
        break;
      }
    }

    if (!is_pivot) continue;

    if (fmpz_is_zero(fmpz_mat_entry(B.mat, row, current_col))) continue;

    pivots.push_back(row);
    fmpz_set(pivot_coeffs + current_col, fmpz_mat_entry(B.mat, row, current_col));
    current_col++;

    if (current_col == b_size) break;
  }

  // Construct matrix of the linear map f acting on B
  fmpq_mat_t f_matrix;
  fmpq_mat_init(f_matrix, b_size, b_size);

  // XXX: in this case, it seems like recomputing f is cheaper than passing in a dense matrix.
  fmpq_mat_t map_of_basis;
  fmpq_mat_init(map_of_basis, b_size, N_basis.size());
  for (int col = 0; col < N_basis.size(); col++) {
    ManinElement fb = f(N_basis[col]);

    auto it = fb.components.begin();
    int pivot_index = 0;

    while (true) {
      if (it == fb.components.end()) break;
      if (pivot_index == pivots.size()) break;

      int row = it->basis_index;

      if (row > pivots[pivot_index]) {
        pivot_index++;
      } else if (row == pivots[pivot_index]) {
        fmpq_t coeff;
        fmpq_init(coeff);
        fmpq_set(coeff, it->coeff);
        fmpq_div_fmpz(coeff, coeff, (pivot_coeffs + pivot_index));
        fmpq_set(fmpq_mat_entry(map_of_basis, pivot_index, col), coeff);
        fmpq_clear(coeff);
        it++;
        pivot_index++;
      } else if (row < pivots[pivot_index]) {
        it++;
      }
    }
  }

  DEBUG_INFO(4,
    {
      printf("map_of_basis: ");
      fmpq_mat_print_dimensions(map_of_basis);
      printf("\n");
    }
  )

  fmpq_mat_mul_fmpz_mat(f_matrix, map_of_basis, B.mat);
  fmpq_mat_clear(map_of_basis);
  flint_cleanup();

  DEBUG_INFO(4,
    {
      printf("f_matrix: ");
      fmpq_mat_print_dimensions(f_matrix);
      printf("\n");
    }
  )

  _fmpz_vec_clear(pivot_coeffs, b_size);

  // Case 1: Space doesn't split: only +1 space
  if (fmpq_mat_is_one(f_matrix)) {
    DEBUG_INFO_PRINT(3, " space doesn't split: A-L action is +1\n");

    fmpq_mat_clear(f_matrix);
    flint_cleanup();

    return {.pos_space = B, .neg_space = dense_empty(N)};
  }

  fmpq_mat_t neg_one;
  fmpq_mat_init(neg_one, b_size, b_size);
  for (int i = 0; i < b_size; i++) {
    fmpq_set_si(fmpq_mat_entry(neg_one, i, i), -1, 1);
  }

  // Case 2: Space doesn't split: only -1 space
  if (fmpq_mat_equal(f_matrix, neg_one)) {
    DEBUG_INFO_PRINT(3, " space doesn't split: A-L action is -1\n");

    fmpq_mat_clear(f_matrix);
    fmpq_mat_clear(neg_one);
    flint_cleanup();

    return {.pos_space = dense_empty(N), .neg_space = B};
  }

  // Case 3: Space splits into +1 and -1 space
  fmpz_mat_t f_matrix_z, kernel, kernel_window, kernel_in_orig_basis;
  fmpz_mat_init(f_matrix_z, b_size, b_size);
  fmpz_mat_init(kernel, b_size, b_size);

  // Set f_matrix to T-1, i.e. positive space
  for (int i = 0; i < b_size; i++) {
    fmpq_sub_si(fmpq_mat_entry(f_matrix, i, i), fmpq_mat_entry(f_matrix, i, i), 1);
  }
  fmpq_mat_get_fmpz_mat_rowwise(f_matrix_z, NULL, f_matrix);

  int rank = fmpz_mat_nullspace_mul(kernel, f_matrix_z);
  fmpz_mat_window_init(kernel_window, kernel, 0, 0, b_size, rank);

  DEBUG_INFO(4,
    {
      printf("(+) kernel_window: ");
      fmpz_mat_print_dimensions(kernel_window);
      printf("\n");
    }
  )

  fmpz_mat_div_colwise_gcd(kernel_window);

  DEBUG_INFO(4,
    {
      printf("(+) kernel_window (cleared): ");
      fmpz_mat_print_dimensions(kernel_window);
      printf("\n");
    }
  )

  fmpz_mat_init(kernel_in_orig_basis, N_basis.size(), rank);
  fmpz_mat_mul(kernel_in_orig_basis, B.mat, kernel_window);
  fmpz_mat_window_clear(kernel_window);
  fmpz_mat_clear(kernel);
  flint_cleanup();

  DEBUG_INFO(4,
    {
      printf("(+) kernel_in_orig_basis: ");
      fmpz_mat_print_dimensions(kernel_in_orig_basis);
      printf("\n");
    }
  )

  fmpz_mat_div_colwise_gcd(kernel_in_orig_basis);

  DEBUG_INFO(4,
    {
      printf("(+) kernel_in_orig_basis (cleared): ");
      fmpz_mat_print_dimensions(kernel_in_orig_basis);
      printf("\n");
    }
  )

  DenseBasis pos_space;
  pos_space.set_move(kernel_in_orig_basis);

  DEBUG_INFO_PRINT(3, " pos_space dimension: %d\n", rank);

  // Set f_matrix from T-1 to T+1, i.e. negative space
  for (int i = 0; i < b_size; i++) {
    fmpq_add_si(fmpq_mat_entry(f_matrix, i, i), fmpq_mat_entry(f_matrix, i, i), 2);
  }

  fmpq_mat_get_fmpz_mat_rowwise(f_matrix_z, NULL, f_matrix);

  fmpz_mat_init(kernel, b_size, b_size);
  rank = fmpz_mat_nullspace_mul(kernel, f_matrix_z);
  fmpz_mat_window_init(kernel_window, kernel, 0, 0, b_size, rank);

  DEBUG_INFO(4,
    {
      printf("(-) kernel_window: ");
      fmpz_mat_print_dimensions(kernel_window);
      printf("\n");
    }
  )

  fmpz_mat_div_colwise_gcd(kernel_window);

  DEBUG_INFO(4,
    {
      printf("(-) kernel_window (cleared): ");
      fmpz_mat_print_dimensions(kernel_window);
      printf("\n");
    }
  )

  fmpz_mat_init(kernel_in_orig_basis, N_basis.size(), rank);
  fmpz_mat_mul(kernel_in_orig_basis, B.mat, kernel_window);
  // fmpz_mat_clear(B_matrix_z);
  fmpz_mat_window_clear(kernel_window);
  fmpz_mat_clear(kernel);
  flint_cleanup();

  DEBUG_INFO(4,
    {
      printf("(-) kernel_in_orig_basis: ");
      fmpz_mat_print_dimensions(kernel_in_orig_basis);
      printf("\n");
    }
  )

  fmpz_mat_div_colwise_gcd(kernel_in_orig_basis);

  DEBUG_INFO(4,
    {
      printf("(-) kernel_in_orig_basis (cleared): ");
      fmpz_mat_print_dimensions(kernel_in_orig_basis);
      printf("\n");
    }
  )

  DenseBasis neg_space;
  neg_space.set_move(kernel_in_orig_basis);

  DEBUG_INFO_PRINT(3, " neg_space dimension: %d\n", rank);

  fmpq_mat_clear(neg_one);
  fmpq_mat_clear(f_matrix);

  fmpz_mat_clear(f_matrix_z);
  // fmpz_mat_clear(kernel);
  fmpz_mat_clear(kernel_in_orig_basis);

  DEBUG_INFO_PRINT(2, "dim %zu -> +: %zu, -: %zu\n", B.mat->c, pos_space.mat->c, neg_space.mat->c);

  return {.pos_space = pos_space, .neg_space = neg_space};
}

DecomposeResult DecomposeResult::empty() {
  return {
    .done = std::vector<Subspace>(),
    .special = std::vector<Subspace>(),
    .remaining = std::vector<Subspace>()
  };
}

// TODO: should be faster to decompose using a matrix of the action on the newform subspace

DecomposeResult decompose(Subspace subspace, FmpqMatrix& map_of_basis, bool dimension_only, bool prime_opt, const slong mem_threshold) {

  auto& B = subspace.basis;

  // If the input space is trivial, the result is also trivial.
  int b_size = B.mat->c;
  if (b_size == 0) {
    return DecomposeResult::empty();
  }
  // Otherwise, we find the level from the first element.
  int64_t N = subspace.level;
  std::vector<ManinBasisElement> N_basis = manin_basis(N);

  // ---------------------------------------- //
  // Setup: constructing matrices for B and f //
  // ---------------------------------------- //

  // Finds pivot rows of the matrix B.
  std::vector<int> pivots;
  fmpz* pivot_coeffs = _fmpz_vec_init(b_size);

  int current_col = 0;
  for (int row = 0; row < N_basis.size(); row++) {
    bool is_pivot = true;
    for (int col = 0; col < b_size; col++) {
      if (col != current_col && !fmpz_is_zero(fmpz_mat_entry(B.mat, row, col))) {
        is_pivot = false;
        break;
      }
    }

    if (!is_pivot) continue;

    if (fmpz_is_zero(fmpz_mat_entry(B.mat, row, current_col))) continue;

    pivots.push_back(row);
    fmpz_set(pivot_coeffs + current_col, fmpz_mat_entry(B.mat, row, current_col));
    current_col++;

    if (current_col == b_size) break;
  }

  // Construct matrix of the linear map f acting on B
  fmpq_mat_t f_matrix;
  fmpq_mat_init(f_matrix, b_size, b_size);

  fmpq_mat_t pivot_rows;
  fmpq_mat_init(pivot_rows, b_size, N_basis.size());

  for (int row = 0; row < b_size; row++) {
    for (int col = 0; col < N_basis.size(); col++) {
      fmpq_div_fmpz(
        fmpq_mat_entry(pivot_rows, row, col),
        fmpq_mat_entry(map_of_basis.mat, pivots[row], col),
        (pivot_coeffs + row)
      );
    }
  }

  // clear map_of_basis after this?

  DEBUG_INFO(4,
    {
      printf("pivot_rows: ");
      fmpq_mat_print_dimensions(pivot_rows);
      printf("\n");
    }
  )

  fmpq_mat_mul_fmpz_mat(f_matrix, pivot_rows, B.mat);
  fmpq_mat_clear(pivot_rows);

  DEBUG_INFO(4,
    {
      printf("f_matrix: ");
      fmpq_mat_print_dimensions(f_matrix);
      printf("\n");
    }
  )

  _fmpz_vec_clear(pivot_coeffs, b_size);

  std::vector<Subspace> done;
  std::vector<Subspace> special;
  std::vector<Subspace> remaining;

  // --------------------------------------------------- //
  // Checks if space doesn't split mod some small primes //
  // --------------------------------------------------- //

  // For spaces of size <= 20, we output the minimal polynomial, so we need to compute it anyways.

  // For spaces of size > 20, we check the factorization of the minimal polynomial mod p, for some
  // small primes p. Any possible factor of the minimal polynomial has degree equal to the sum of
  // the degrees of some factors mod p, so by checking various p we can potentially deduce that
  // the minimal polynomial is irreducible.
  if (b_size > 20 && prime_opt) {
    fmpz_mat_t f_matrix_z;
    fmpz_mat_init(f_matrix_z, b_size, b_size);

    fmpz_t temp_den;
    fmpz_init(temp_den);
    fmpq_mat_get_fmpz_mat_matwise(f_matrix_z, temp_den, f_matrix);
    fmpz_clear(temp_den);

    DEBUG_INFO(5,
      {
        printf("f_matrix_z: ");
        fmpz_mat_print_dimensions(f_matrix_z);
        printf("\n");
      }
    )

    ulong p = 32;
    std::vector<bool> possible_factor_degs(b_size + 1);
    for (int i = 0; i <= b_size; i++)
      possible_factor_degs[i] = true;

    for (int iter = 0; iter < 10; iter++) {
      p = n_nextprime(p, 1);

      nmod_mat_t f_matrix_mod_p;
      nmod_mat_init(f_matrix_mod_p, b_size, b_size, p);

      nmod_poly_t min_poly_mod_p;
      nmod_poly_init(min_poly_mod_p, p);

      fmpz_mat_get_nmod_mat(f_matrix_mod_p, f_matrix_z);

      DEBUG_INFO_PRINT(5, "f_matrix_mod_p\n");

      nmod_mat_minpoly(min_poly_mod_p, f_matrix_mod_p);

      DEBUG_INFO_PRINT(5, "min_poly_mod_p deg: %d\n", min_poly_mod_p->length - 1);

      if (min_poly_mod_p->length == b_size + 1) {
        // XXX: Minimal polynomial is rarely irreducible, so we don't save time by checking irreducibility first.
        //
        // if (nmod_poly_is_irreducible(min_poly_mod_p)) {
        //   done.push_back(B);
        //   DEBUG_INFO_PRINT(3, " minimal polynomial irreducible and degree equal to space dimension %zu in mod %ld\n", b_size, p);

        //   nmod_mat_clear(f_matrix_mod_p);
        //   nmod_poly_clear(min_poly_mod_p);

        //   fmpq_mat_clear(f_matrix);
        //   fmpz_mat_clear(B.mat);
        //   return {.done = done, .special = special, .remaining = remaining};
        // } else {
        nmod_poly_factor_t facs;
        nmod_poly_factor_init(facs);
        nmod_poly_factor(facs, min_poly_mod_p);

        DEBUG_INFO(4,
          {
            printf("p = %ld, deg = %ld, facs = ", p, min_poly_mod_p->length - 1);
            for (int i = 0; i < facs->num; i++) {
              printf("%ld,", (facs->p + i)->length - 1);
            }
            printf("\n");
          }
        )

        std::vector<bool> p_factor_degs(b_size + 1);
        p_factor_degs[0] = true;
        for (int i = 1; i <= b_size; i++)
          p_factor_degs[i] = false;

        int sum = 0;
        for (int i = 1; i <= facs->num; i++) {
          int factor_deg = (facs->p + (i % facs->num))->length - 1;
          sum += factor_deg;
          for (int j = sum; j >= factor_deg; j--) {
            p_factor_degs[j] = p_factor_degs[j] || p_factor_degs[j - factor_deg];
          }
        }

        bool irred = true;
        for (int j = 1; j < b_size; j++) {
          possible_factor_degs[j] = possible_factor_degs[j] && p_factor_degs[j];
          irred = irred && !(possible_factor_degs[j]);
        }

        nmod_poly_factor_clear(facs);

        if (irred) {
          done.push_back(subspace);

          DEBUG_INFO_PRINT(3, " minimal polynomial irreducible and degree equal to space dimension %zu, found at iteration %ld\n", b_size, iter);

          nmod_mat_clear(f_matrix_mod_p);
          nmod_poly_clear(min_poly_mod_p);

          fmpq_mat_clear(f_matrix);
          // fmpz_mat_clear(B_matrix_z);
          fmpz_mat_clear(f_matrix_z);
          return {.done = done, .special = special, .remaining = remaining};
        }
      }
      else {
        DEBUG_INFO_PRINT(4, "p = %ld, deg = %ld\n", p, min_poly_mod_p->length - 1);
      }


      nmod_mat_clear(f_matrix_mod_p);
      nmod_poly_clear(min_poly_mod_p);
    }

    fmpz_mat_clear(f_matrix_z);
  }

  // --------------------- //
  // Computes minpoly of f //
  // --------------------- //

  // Note: computing char_poly seems much slower than min_poly.
  fmpq_poly_t min_poly;
  fmpz_poly_t min_poly_z;
  fmpq_poly_init(min_poly);
  fmpz_poly_init(min_poly_z);

  if (fmpq_mat_is_zero(f_matrix)) {
    // XXX: FLINT seems to think that the minimal polynomial of the zero matrix is 1, and not T.
    // Manually set the minimal polynomial to x.
    fmpz_poly_set_coeff_si(min_poly_z, 1, 1);
  } else {
    fmpq_mat_minpoly(min_poly, f_matrix);
    fmpq_poly_get_numerator(min_poly_z, min_poly);
  }

  int min_poly_degree = fmpz_poly_degree(min_poly_z);

  DEBUG_INFO_PRINT(3, " min_poly_z degree: %d\n", min_poly_degree);

  DEBUG_INFO(5,
    {
      printf(" min_poly_z: ");
      fmpz_poly_print_pretty(min_poly_z, "T");
      printf("\n");
    }
  )

  fmpq_poly_clear(min_poly);

  // ---------------------------- //
  // Computes space decomposition //
  // ---------------------------- //

  fmpz_poly_factor_t min_poly_factored;
  fmpz_poly_factor_init(min_poly_factored);
  fmpz_poly_factor(min_poly_factored, min_poly_z);

  int num_factors = min_poly_factored->num;
  if (num_factors == 0) {
    // This should actually be impossible
    assert(false);
  } else if (num_factors == 1) {
    int deg = fmpz_poly_degree(min_poly_factored->p);
    if (deg == b_size) {
      Subspace out = subspace;
      FmpzPoly out_poly;
      out_poly.set_copy(min_poly_z);
      out.hecke_field_poly = out_poly;
      done.emplace_back(out);
      DEBUG_INFO_PRINT(3, " minimal polynomial irreducible and degree equal to space dimension %zu\n", b_size);
    } else {
      special.push_back(subspace);
      DEBUG_INFO_PRINT(3, " minimal polynomial irreducible but space dimension %zu is not equal to degree %d\n", b_size, deg);
    }
  } else {
    fmpq_mat_t poly_on_f_matrix;
    fmpz_mat_t poly_on_f_matrix_z, poly_mat_kernel;
    fmpq_mat_init(poly_on_f_matrix, b_size, b_size);
    fmpz_mat_init(poly_on_f_matrix_z, b_size, b_size);
    fmpz_mat_init(poly_mat_kernel, b_size, b_size);

    fmpz_mat_t poly_mat_kernel_window, poly_mat_kernel_in_orig_basis;

    int dimension_excess = 0;

    // For each factor g of the minpoly, we compute g(T) and take its kernel.
    for (int i = 0; i < num_factors; i++) {

      fmpz_poly_struct *factor = min_poly_factored->p + i;
      int exp = *(min_poly_factored->exp + i);
      int degree = fmpz_poly_degree(factor);

      if (dimension_only && (min_poly_degree + dimension_excess + degree > b_size)) {
        DEBUG_INFO_PRINT(3, " factor appears only once, degree: %d\n", degree);
        fmpz_mat_t out_mat;
        fmpz_mat_init(out_mat, N_basis.size(), degree);
        DenseBasis output;
        output.set_move(out_mat);
        done.emplace_back(output, true, N, subspace.atkin_lehner_pos, subspace.atkin_lehner_neg);
        continue;
      }

      fmpz_poly_apply_fmpq_mat(poly_on_f_matrix, f_matrix, factor, mem_threshold);

      DEBUG_INFO(4,
        {
          printf("poly_on_f_matrix: ");
          fmpq_mat_print_dimensions(poly_on_f_matrix);
          printf("\n");
        }
      )
      // NOTE: this needs to be rowwise!
      fmpq_mat_get_fmpz_mat_rowwise(poly_on_f_matrix_z, NULL, poly_on_f_matrix);

      DEBUG_INFO(4,
        {
          printf("poly_on_f_matrix_z: ");
          fmpz_mat_print_dimensions(poly_on_f_matrix_z);
          printf("\n");
        }
      )

      int rank = fmpz_mat_nullspace_mul(poly_mat_kernel, poly_on_f_matrix_z);

      fmpz_mat_window_init(poly_mat_kernel_window, poly_mat_kernel, 0, 0, b_size, rank);

      DEBUG_INFO(4,
        {
          printf("poly_mat_kernel_window: ");
          fmpz_mat_print_dimensions(poly_mat_kernel_window);
          printf("\n");
        }
      )

      fmpz_mat_div_colwise_gcd(poly_mat_kernel_window);

      DEBUG_INFO(6,
        {
          fmpz_mat_print_pretty(poly_mat_kernel_window);
          printf("\n");
        }
      )

      DEBUG_INFO(4,
        {
          printf("poly_mat_kernel_window (cleared): ");
          fmpz_mat_print_dimensions(poly_mat_kernel_window);
          printf("\n");
        }
      )

      fmpz_mat_init(poly_mat_kernel_in_orig_basis, N_basis.size(), rank);
      fmpz_mat_mul(poly_mat_kernel_in_orig_basis, B.mat, poly_mat_kernel_window);
      fmpz_mat_window_clear(poly_mat_kernel_window);

      DEBUG_INFO(4,
        {
          printf("poly_mat_kernel_in_orig_basis: ");
          fmpz_mat_print_dimensions(poly_mat_kernel_in_orig_basis);
          printf("\n");
        }
      )

      fmpz_mat_div_colwise_gcd(poly_mat_kernel_in_orig_basis);

      DEBUG_INFO(4,
        {
          printf("poly_mat_kernel_in_orig_basis (cleared): ");
          fmpz_mat_print_dimensions(poly_mat_kernel_in_orig_basis);
          printf("\n");
        }
      )

      DenseBasis output;
      output.set_move(poly_mat_kernel_in_orig_basis);

      DEBUG_INFO(3,
        {
          printf(" subspace dimension: %d, factor degree: %d\n", rank, degree);
        }
      )

      if (degree == rank) {
        Subspace out(output, true, N, subspace.atkin_lehner_pos, subspace.atkin_lehner_neg);
        FmpzPoly out_poly;
        out_poly.set_copy(factor);
        out.hecke_field_poly = out_poly;
        done.emplace_back(out);
      } else {
        dimension_excess += (rank - degree);
        Subspace out(output, false, N, subspace.atkin_lehner_pos, subspace.atkin_lehner_neg);
        remaining.emplace_back(out);
      }
    }

    fmpq_mat_clear(poly_on_f_matrix);
    fmpz_mat_clear(poly_on_f_matrix_z);
    fmpz_mat_clear(poly_mat_kernel);
  }

  fmpz_poly_clear(min_poly_z);
  fmpz_poly_factor_clear(min_poly_factored);

  fmpq_mat_clear(f_matrix);
  // fmpz_mat_clear(B_matrix_z);

  return {.done = done, .special = special, .remaining = remaining};
}

