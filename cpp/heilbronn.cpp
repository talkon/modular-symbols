#include "heilbronn.h"
#include <cmath>

// Algorithm given in Cremona Ch 2.4
std::vector<IntMatrix2x2> heilbronn_cremona(int64_t p) {
  std::vector<IntMatrix2x2> result;
  result.push_back({.x = 1, .y = 0, .z = 0, .w = p});

  for (int s = 0; s < p; s++) {
    int r = s - (p - 1) / 2;
    int x1 = p;
    int x2 = -r;
    int y1 = 0;
    int y2 = 1;
    int a = -p;
    int b = r;
    result.push_back({.x = x1, .y = x2, .z = y1, .w = y2});
    while (b != 0) {
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

// Merel's set X (see p.87 in Merel)
// Based on algorithm used in Sage at src/sage/modular/modsym/heilbronn.pyx in the Sage source code.
std::vector<IntMatrix2x2> heilbronn_merel(int64_t n) {
  std::vector<IntMatrix2x2> result;

  int a, b, c, d, q;
  for (a = 1; a <= n; a++) {
    q = n / a;
    if (q * a == n) {
      d = q;
      for (b = 0; b < a; b++) {
        result.push_back({.x = a, .y = b, .z = 0, .w = d});
      }
      for (c = 1; c < d; c++) {
        result.push_back({.x = a, .y = 0, .z = c, .w = d});
      }
    }
    for (d = q + 1; d <= n; d++) {
      int bc = a * d - n;
      for (c = bc / a + 1; c < d; c++) {
        if (bc % c == 0) {
          result.push_back({.x = a, .y = bc/c, .z = c, .w = d});
        }
      }
    }
  }

  return result;
}