#include "../cpp/newspace.h"

#include <iostream>
#include <gtest/gtest.h>

void expect_newspace_basis_size(int level, int size) {
  std::vector<ManinElement> basis = newspace_basis(level);
  EXPECT_EQ(basis.size(), size);
}

TEST(NewspaceTests, BasisSizes_Small) {
  expect_newspace_basis_size(95, 7);
  expect_newspace_basis_size(96, 2);
  expect_newspace_basis_size(97, 7);
  expect_newspace_basis_size(98, 3);
  expect_newspace_basis_size(99, 4);
  expect_newspace_basis_size(100, 1);
}

TEST(NewspaceTests, BasisSizes_Composite) {
  expect_newspace_basis_size(400, 8);
  expect_newspace_basis_size(420, 4);
  expect_newspace_basis_size(462, 9);
}

TEST(NewspaceTests, BasisSizes_Semiprime) {
  expect_newspace_basis_size(403, 31);
  expect_newspace_basis_size(437, 33);
  expect_newspace_basis_size(451, 33);
}

TEST(NewspaceTests, BasisSizes_Prime) {
  expect_newspace_basis_size(401, 33);
  expect_newspace_basis_size(421, 34);
  expect_newspace_basis_size(463, 38);
}

TEST(NewspaceTests, BasisSizes_2532) {
  expect_newspace_basis_size(2532, 36);
}