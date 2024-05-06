#include "modular_symbol.h"
#include "manin_symbol.h"
#include "manin_basis.h"
#include "manin_element.h"
#include "boundary_map.h"
#include "newspace.h"
#include "newform_subspaces.h"
#include "debug_utils.h"

#include <unistd.h>

void print_help() {
  printf("Usage\
          \n -n N : (required) sets level to N\
          \n -t T : sets trace depth to T (default 10)\
          \n -s   : compute only subspace dimensions (overrides -t)\
          \n -p   : disable newest optimization (for testing only)\
          \n -v V : sets verbosity to V\
          \n -d D : (used for debugging purposes)\
          \n -h   : prints this help message and exits\
          \n");
}

int main(int argc, char** argv) {

  char c;

  int64_t level = 0;
  int verbose = 0;
  int trace_depth = 10;
  bool dimension_only = false;
  bool prime_opt = true;

  if (argc == 1) {
    print_help();
    return 0;
  }

  while ((c = getopt (argc, argv, "n:v:d:t:hsp")) != -1) {
    switch (c) {
      case 'n':
        level = atol(optarg);
        break;
      case 'v':
        verbose = atoi(optarg);
        break;
      case 't':
        trace_depth = atoi(optarg);
        break;
      case 's':
        dimension_only = true;
        break;
      case 'd':
        set_verbosity(10);
        debug_temp(atoi(optarg));
        return 0;
      case 'p':
        prime_opt = false;
        break;
      case 'h':
        print_help();
        return 0;
    }
  }

  if (dimension_only) trace_depth = 0;

  set_verbosity(verbose);
  init_timer();

  DEBUG_INFO_PRINT(1, "Started computation for level %lld\n", level);
  auto subspaces = newform_subspaces(level, dimension_only, trace_depth, prime_opt);
  DEBUG_INFO_PRINT(1, "Finished computation for level %lld\n", level);

  for (auto& subspace : subspaces) {
    subspace.print();
    printf("\n");
  }

  printf("%lld:[", level);
  for (int i = 0; i < subspaces.size(); i++) {
    printf("%d", subspaces[i].dimension());
    if (i != subspaces.size() - 1) printf(",");
  }
  printf("]:%.4lf\n", get_elapsed_time());
  return 0;
}