#ifndef DEBUG_UTILS_H
#define DEBUG_UTILS_H

#include <iostream>
#include <chrono>

// --- Printing ---

// color macros, used to color debugging messages
#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

// macro to probe how far a program was executed before breaking
#define DEBUG_PRINT printf("function " MAG"%s() " RESET "| line " CYN "%d\n" RESET, __func__, __LINE__);

// --- Timing ---

inline std::chrono::system_clock::time_point start_time;

inline void init_timer() {
  start_time = std::chrono::system_clock::now();
}

inline double get_elapsed_time() {
  std::chrono::duration<double> elapsed = (std::chrono::system_clock::now() - start_time);
  return elapsed.count();
}

inline void info_with_time() {
  printf(BLU "%.4lfs > " RESET, get_elapsed_time());
}

#define DEBUG_INFO_PRINT(v, ...) if (verbosity >= v) { info_with_time(); printf(__VA_ARGS__); }

#define DEBUG_INFO(v, ...) if (verbosity >= v) { info_with_time(); __VA_ARGS__ }

// --- Verbosity level ---

inline int verbosity = 0;

inline void set_verbosity(int v) {
  verbosity = v;
}

// --- Functions used for debugging ---

// Debugging function that is run when using `./bin/main -d`.
void debug_temp(int d);

// Attempts a 4x3 matrix nullspace computation.
// Can detect issues when FLINT is trying to access an mpz that has already been freed.
void check_status();

// Probes FLINT's mpz freelist by requesting `depth` new mpzs.
// Can detect double-frees (which will end up duplicated on the free list).
void probe_fmpz_freelist(int depth);

#endif // DEBUG_UTILS_H