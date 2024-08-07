#include "modular_symbol.h"
#include "manin_symbol.h"
#include "manin_basis.h"
#include "manin_element.h"
#include "boundary_map.h"
#include "newspace.h"
#include "newform_subspaces.h"
#include "debug_utils.h"

#include <flint/flint.h>
#include <unistd.h>

void print_help() {
  printf("Usage\
          \n -n <level>           : (required) sets level to <level>\
          \n -t <min_trace_depth> : sets minimum trace depth (default 0), i.e. \
          \n                      : will always compute the trace form to at least this depth\
          \n -T <max_trace_depth> : sets maximum trace depth (default unbounded), i.e.\
          \n                        will not compute beyond this depth even if some subspaces\
          \n                        still have the same trace form up to <max_trace_depth>\
          \n -M <mem_threshold>   : (experimental, only on Linux) sets memory threshold in matrix\
          \n                        polynomial evaluation, in MB. Empirically, use value of m to keep\
          \n                        overall memory usage below 2.5m _if_ this step is the bottleneck\
          \n -s                   : compute only subspace dimensions (overrides -t and -T)\
          \n -v <verbosity>       : sets verbosity to V\
          \n -p                   : disable newest optimization (used for debugging)\
          \n -d <debug_arg>       : (used for debugging)\
          \n -h                   : prints this help message and exits\
          \n");
}

int main(int argc, char** argv) {

  char c;

  int64_t level = 0;
  int verbose = 0;
  int min_trace_depth = 0;
  int max_trace_depth = -1;
  slong mem_threshold = 0;
  bool dimension_only = false;
  bool prime_opt = true;

  if (argc == 1) {
    print_help();
    return 0;
  }

  while ((c = getopt (argc, argv, "n:v:d:t:T:M:hsp")) != -1) {
    switch (c) {
      case 'n':
        level = atol(optarg);
        break;
      case 'v':
        verbose = atoi(optarg);
        break;
      case 't':
        min_trace_depth = atoi(optarg);
        break;
      case 'T':
        max_trace_depth = atoi(optarg);
        break;
      case 'M':
        mem_threshold = atol(optarg) << 10L;
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

  if (dimension_only) {
    min_trace_depth = 0;
    max_trace_depth = 0;
  }

  set_verbosity(verbose);
  init_timer();

  DEBUG_INFO_PRINT(1, "Started computation for level %lld\n", level);
  auto subspaces = newform_subspaces(level, dimension_only, min_trace_depth, max_trace_depth, prime_opt, mem_threshold);
  DEBUG_INFO_PRINT(1, "Finished computation for level %lld\n", level);

  int si = 0;
  for (auto& subspace : subspaces) {
    subspace.print(si);
    printf("\n");
    si++;
  }

  printf("%lld:[", level);
  for (int i = 0; i < subspaces.size(); i++) {
    printf("%d", subspaces[i].dimension());
    if (i != subspaces.size() - 1) printf(",");
  }
  printf("]:%.4lf\n", get_elapsed_time());
  flint_cleanup_master();
  return 0;
}
