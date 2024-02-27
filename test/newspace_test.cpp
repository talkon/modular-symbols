#include "../cpp/newspace.h"

#include <iostream>
#include <gtest/gtest.h>

void expect_newspace_basis_size(int level, int size) {
  std::vector<ManinElement> basis = newspace_basis(level);
  EXPECT_EQ(basis.size(), size);
}

TEST(NewspaceTests, BasisSizes) {
  expect_newspace_basis_size(95, 7);
  expect_newspace_basis_size(96, 2);
  expect_newspace_basis_size(97, 7);
  expect_newspace_basis_size(98, 3);
  expect_newspace_basis_size(99, 4);
  expect_newspace_basis_size(100, 1);

  expect_newspace_basis_size(400, 8);
  expect_newspace_basis_size(420, 4);
  expect_newspace_basis_size(462, 9);
}