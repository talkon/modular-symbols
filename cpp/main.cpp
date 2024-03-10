#include "modular_symbol.h"
#include "manin_symbol.h"
#include "manin_basis.h"
#include "manin_element.h"
#include "boundary_map.h"
#include "newspace.h"
#include "newform_subspaces.h"
#include "debug_utils.h"

#include <unistd.h>

int main(int argc, char** argv) {

  char c;

  int64_t level = 0;
  int verbose = 0;
  bool use_atkin_lehner = false;

  while ((c = getopt (argc, argv, "n:v:ah")) != -1) {
    switch (c) {
      case 'n':
        level = atol(optarg);
        break;
      case 'v':
        verbose = atoi(optarg);
        break;
      case 'a':
        use_atkin_lehner = true;
        break;
      case 'h':
        printf("Usage\
        \n -n N : (required) sets level to N\
        \n -v V : sets verbosity to V\
        \n -a   : uses Atkin-Lehner involutions\
        \n -h   : prints this help message and exits\
        \n");
        return 0;
    }
  }

  set_verbosity(verbose);
  init_timer();

  auto dims = newform_subspace_dimensions(level, use_atkin_lehner);
  DEBUG_INFO(1, "finished computation, output:\n");
  printf("%lld:[", level);
  for (int i = 0; i < dims.size(); i++) {
    printf("%d", dims[i]);
    if (i != dims.size() - 1) printf(",");
  }
  printf("]:%.4lf\n", get_elapsed_time());
  return 0;
}