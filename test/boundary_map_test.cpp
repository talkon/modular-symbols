#include "../cpp/boundary_map.h"

#include <iostream>
#include <gtest/gtest.h>

void expect_cuspidal_basis_size(int level, int size) {
  std::vector<ManinElement> basis = cuspidal_manin_basis(level);
  EXPECT_EQ(basis.size(), size);
}

TEST(BoundaryMapTest, BasisSizes_Small) {
  expect_cuspidal_basis_size(95, 9);
  expect_cuspidal_basis_size(96, 9);
  expect_cuspidal_basis_size(97, 7);
  expect_cuspidal_basis_size(98, 7);
  expect_cuspidal_basis_size(99, 9);
  expect_cuspidal_basis_size(100, 7);
}

TEST(BoundaryMapTest, BasisSizes_Composite) {
  expect_cuspidal_basis_size(400, 43);
  expect_cuspidal_basis_size(420, 85);
  expect_cuspidal_basis_size(462, 89);
}

TEST(BoundaryMapTest, BasisSizes_Semiprime) {
  expect_cuspidal_basis_size(403, 35);
  expect_cuspidal_basis_size(437, 39);
  expect_cuspidal_basis_size(451, 41);
}

TEST(BoundaryMapTest, BasisSizes_Prime) {
  expect_cuspidal_basis_size(401, 33);
  expect_cuspidal_basis_size(421, 34);
  expect_cuspidal_basis_size(463, 38);
}