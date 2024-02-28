#include "../cpp/newform_subspaces.h"

#include <iostream>
#include <gtest/gtest.h>

void expect_newform_subspace_dimensions(int level, std::vector<int> expected_sizes) {
  std::vector<int> computed_sizes = newform_subspace_dimensions(level);
  EXPECT_EQ(computed_sizes, expected_sizes);
}

TEST(NewformSubspaceTests, BasisSizes_Small) {
  expect_newform_subspace_dimensions(95, {3, 4});
  expect_newform_subspace_dimensions(96, {1, 1});
  expect_newform_subspace_dimensions(97, {3, 4});
  expect_newform_subspace_dimensions(98, {1, 2});
  expect_newform_subspace_dimensions(99, {1, 1, 1, 1});
  expect_newform_subspace_dimensions(100, {1});
}

TEST(NewformSubspaceTests, BasisSizes_Composite) {
  expect_newform_subspace_dimensions(400, {1, 1, 1, 1, 1, 1, 1, 1});
  expect_newform_subspace_dimensions(420, {1, 1, 1, 1});
  expect_newform_subspace_dimensions(462, {1, 1, 1, 1, 1, 1, 1, 2});
}

TEST(NewformSubspaceTests, BasisSizes_Semiprime) {
  expect_newform_subspace_dimensions(403, {2, 6, 7, 8, 8});
  expect_newform_subspace_dimensions(437, {1, 1, 2, 2, 2, 5, 8, 12});
  expect_newform_subspace_dimensions(451, {1, 5, 5, 10, 12});
}

TEST(NewformSubspaceTests, BasisSizes_Prime) {
  expect_newform_subspace_dimensions(401, {12, 21});
  expect_newform_subspace_dimensions(421, {15, 19});
  expect_newform_subspace_dimensions(463, {16, 22});
}