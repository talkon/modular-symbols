#include "manin_symbol.h"
#include "manin_basis.h"
#include "boundary_map.h"

int main(int argc, char** argv) {
  // // Tests manin_generators
  // std::vector<ManinGenerator> mgs;
  // for (int i = 0; i < 10; i++) {
  //   mgs = manin_generators(2 * 3 * 5 * 7 * 11 * 13);
  //   printf("%zu\n", mgs.size());
  // }
  // mgs[2 * 3 * 5 * 7 * 11 * 13 + 50000].print();
  // printf("\n");

  // // Tests find_generator
  // for (int i = 11; i < 20; i++) {
  //   ManinSymbol ms1 = {.N = i, .c = 4, .d = 7};
  //   ManinGenerator mg1 = find_generator(ms1);
  //   ms1.print();
  //   mg1.print();
  //   printf(" %lld %p\n", mg1.index, &mg1);
  // }

  // // Tests relation matrix
  // int level = atoi(argv[1]);
  // std::vector<ManinGenerator> basis = manin_basis(level);

  // printf("[output] manin_basis size: %zu, basis:\n", basis.size());
  // for (ManinGenerator generator : basis) {
  //   generator.print();
  //   printf("\n");
  // }

  // // Tests boundary map
  // std::vector<ManinElement> cuspidal_basis = cuspidal_manin_basis(level);
  // printf("[output] cuspidal_basis size: %zu, basis:\n", cuspidal_basis.size());
  // for (ManinElement element : cuspidal_basis) {
  //   element.print_with_generators();
  //   printf("\n");
  // }

  // Tests continued fractions
  fraction_to_manin_element(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));

  flint_cleanup_master();

  return 0;
}