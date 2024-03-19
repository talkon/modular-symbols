#ifndef NEWSPACE_H
#define NEWSPACE_H

#include "manin_element.h"

#include <vector>

// Computes a basis of the newspace of a given level.
std::vector<ManinElement> newspace_basis(int64_t);

#endif // NEWSPACE_H
