#include "heilbronn.h"

std::vector<IntMatrix2x2> heilbronn_cremona(int64_t p) {
  std::vector<IntMatrix2x2> result;
  result.push_back({.x = 1, .y = 0, .z = 0, .w = p});

  for (int r = 0; r < p; r++) {
    int x1 = p;
    int x2 = -r;
    int y1 = 0;
    int y2 = 1;
    int a = -p;
    int b = r;
    result.push_back({.x = x1, .y = x2, .z = y1, .w = y2});
    while (b != 0) {
      // XXX: This rounding behavior is slightly inefficient (gives slightly more matrices than necessary), but is correct afaict.
      // TODO: fix this rounding
      int q = std::lround((double) a / (double) b);
      int c = a - b * q;
      a = -b;
      b = c;
      int x3 = q * x2 - x1;
      x1 = x2;
      x2 = x3;
      int y3 = q * y2 - y1;
      y1 = y2;
      y2 = y3;
      result.push_back({.x = x1, .y = x2, .z = y1, .w = y2});
    }
  }

  return result;
}