#include "../cpp/manin_basis.h"

#include <iostream>
#include <gtest/gtest.h>

void expect_basis_size(int level, int size) {
  std::vector<ManinBasisElement> basis = manin_basis(level);
  EXPECT_EQ(basis.size(), size);
}

TEST(ManinBasisTest, BasisSizes_Small) {
  expect_basis_size(95, 12);
  expect_basis_size(96, 20);
  expect_basis_size(97, 8);
  expect_basis_size(98, 16);
  expect_basis_size(99, 14);
  expect_basis_size(100, 18);
}

TEST(ManinBasisTest, BasisSizes_Composite) {
  expect_basis_size(400, 64);
  expect_basis_size(420, 108);
  expect_basis_size(462, 104);
}

TEST(ManinBasisTest, BasisSizes_Semiprime) {
  expect_basis_size(403, 38);
  expect_basis_size(437, 42);
  expect_basis_size(451, 44);
}

TEST(ManinBasisTest, BasisSizes_Prime) {
  expect_basis_size(401, 34);
  expect_basis_size(421, 35);
  expect_basis_size(463, 39);
}

TEST(ManinBasisTest, ResultIsCached) {
  for (int i = 0; i < 10000; i++) {
    expect_basis_size(210, 56);
  }
}
