#ifndef NEWSPACE_H
#define NEWSPACE_H

#include <vector>

// Forward declaration
class ManinElement;

// Computes a basis of the newspace of a given level.
std::vector<ManinElement> newspace_basis(int64_t);

#endif // NEWSPACE_H
