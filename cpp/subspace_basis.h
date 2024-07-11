#ifndef SUBSPACE_BASIS_H
#define SUBSPACE_BASIS_H

#include "manin_element.h"
#include "flint_wrappers.h"

#include <vector>

typedef std::vector<ManinElement> SparseBasis;
typedef FmpzMatrix DenseBasis;

DenseBasis sparse_to_dense(SparseBasis& sparse, int64_t level);

SparseBasis dense_to_sparse(DenseBasis& dense, int64_t level);

#endif // SUBSPACE_BASIS_H