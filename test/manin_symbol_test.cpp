#include "../cpp/manin_symbol.h"

#include <gtest/gtest.h>

TEST(ManinBasis, BasisSizes) {
  std::vector<ManinGenerator> basis = manin_basis(462);
  EXPECT_EQ(basis.size(), 104);
}
