#ifndef LINALG_H
#define LINALG_H

// This file defines linear algebra helper functions

#include "manin_basis.h"
#include "manin_element.h"

#include <vector>

// Computes the kernel (represented by a basis) of a given linear map `f` acting on
// a vector space of ManinElements (also represented by a basis).
// `f` should be a map to Manin symbols of level `M`.
std::vector<ManinElement> map_kernel(std::vector<ManinElement>, std::function<ManinElement(ManinBasisElement)>, int64_t M);

#endif // LINALG_H