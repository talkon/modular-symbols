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

  while ((c = getopt (argc, argv, "n:v:hd")) != -1) {
    switch (c) {
      case 'n':
        level = atol(optarg);
        break;
      case 'v':
        verbose = atoi(optarg);
        break;
      case 'd':
        set_verbosity(10);
        debug_temp();
        return 0;
      case 'h':
        printf("Usage\
        \n -n N : (required) sets level to N\
        \n -v V : sets verbosity to V\
        \n -d   : (used for debugging purposes)\
        \n -h   : prints this help message and exits\
        \n");
        return 0;
    }
  }

  set_verbosity(verbose);
  init_timer();

  DEBUG_INFO_PRINT(1, "Started computation for level %lld\n", level);
  auto dims = newform_subspace_dimensions(level);
  DEBUG_INFO_PRINT(1, "Finished computation for level %lld\n", level);
  printf("%lld:[", level);
  for (int i = 0; i < dims.size(); i++) {
    printf("%d", dims[i]);
    if (i != dims.size() - 1) printf(",");
  }
  printf("]:%.4lf\n", get_elapsed_time());
  return 0;
}