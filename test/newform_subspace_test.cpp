#include "../cpp/newform_subspaces.h"

#include <iostream>
#include <gtest/gtest.h>

void expect_newform_subspace_dimensions(int level, bool use_atkin_lehner, std::vector<int> expected_sizes) {
  std::vector<int> computed_sizes = newform_subspace_dimensions(level, use_atkin_lehner);
  EXPECT_EQ(computed_sizes, expected_sizes);
}

TEST(NewformSubspaceTests, BasisSizes_Small) {
  expect_newform_subspace_dimensions(95, false, {3, 4});
  expect_newform_subspace_dimensions(96, false, {1, 1});
  expect_newform_subspace_dimensions(97, false, {3, 4});
  expect_newform_subspace_dimensions(98, false, {1, 2});
  expect_newform_subspace_dimensions(99, false, {1, 1, 1, 1});
  expect_newform_subspace_dimensions(100, false,  {1});
}

TEST(NewformSubspaceTests, BasisSizes_Composite) {
  expect_newform_subspace_dimensions(400, false, {1, 1, 1, 1, 1, 1, 1, 1});
  expect_newform_subspace_dimensions(420, false, {1, 1, 1, 1});
  expect_newform_subspace_dimensions(462, false, {1, 1, 1, 1, 1, 1, 1, 2});
}

TEST(NewformSubspaceTests, BasisSizes_Composite_AL) {
  expect_newform_subspace_dimensions(400, true, {1, 1, 1, 1, 1, 1, 1, 1});
  expect_newform_subspace_dimensions(420, true, {1, 1, 1, 1});
  expect_newform_subspace_dimensions(462, true, {1, 1, 1, 1, 1, 1, 1, 2});
}

TEST(NewformSubspaceTests, BasisSizes_Semiprime) {
  expect_newform_subspace_dimensions(403, false, {2, 6, 7, 8, 8});
  expect_newform_subspace_dimensions(437, false, {1, 1, 2, 2, 2, 5, 8, 12});
  expect_newform_subspace_dimensions(451, false, {1, 5, 5, 10, 12});
}

TEST(NewformSubspaceTests, BasisSizes_Semiprime_AL) {
  expect_newform_subspace_dimensions(403, true, {2, 6, 7, 8, 8});
  expect_newform_subspace_dimensions(437, true, {1, 1, 2, 2, 2, 5, 8, 12});
  expect_newform_subspace_dimensions(451, true, {1, 5, 5, 10, 12});
}

TEST(NewformSubspaceTests, BasisSizes_Prime) {
  expect_newform_subspace_dimensions(401, false, {12, 21});
  expect_newform_subspace_dimensions(421, false, {15, 19});
  expect_newform_subspace_dimensions(463, false, {16, 22});
}

TEST(NewformSubspaceTests, BasisSizes_Prime_AL) {
  expect_newform_subspace_dimensions(401, true, {12, 21});
  expect_newform_subspace_dimensions(421, true, {15, 19});
  expect_newform_subspace_dimensions(463, true, {16, 22});
}