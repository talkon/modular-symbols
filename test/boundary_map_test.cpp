#include "../cpp/boundary_map.h"

#include <iostream>
#include <gtest/gtest.h>

void expect_cuspidal_basis_size(int level, int size) {
  std::vector<ManinElement> basis = cuspidal_manin_basis(level);
  EXPECT_EQ(basis.size(), size);
}

TEST(BoundaryMapTest, BasisSizes) {
  expect_cuspidal_basis_size(95, 9);
  expect_cuspidal_basis_size(96, 9);
  expect_cuspidal_basis_size(97, 7);
  expect_cuspidal_basis_size(98, 7);
  expect_cuspidal_basis_size(99, 9);
  expect_cuspidal_basis_size(100, 7);

  expect_cuspidal_basis_size(400, 43);
  expect_cuspidal_basis_size(420, 85);
  expect_cuspidal_basis_size(462, 89);
}