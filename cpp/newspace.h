#ifndef NEWSPACE_H
#define NEWSPACE_H

#include <vector>
#include "manin_element.h"

// Computes a basis of the newspace of a given level.
std::vector<ManinElement> newspace_basis(int64_t);

#endif // NEWSPACE_H
