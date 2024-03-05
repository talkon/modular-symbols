#ifndef DEBUG_TIMER_H
#define DEBUG_TIMER_H

#include <iostream>
#include <chrono>

inline std::chrono::system_clock::time_point start_time;

inline void init_debug_time() {
  start_time = std::chrono::system_clock::now();
}

inline void info_with_time() {
  auto now = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = now - start_time;
  printf("[info] %7.4lfs >", elapsed_seconds);
}

#endif // DEBUG_TIMER_H