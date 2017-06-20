#ifndef BSVD_OMP_H
#define BSVD_OMP_H

#include "binmat.h"

void initialize_model(const binary_matrix& E, 
		      binary_matrix& D, 
		      binary_matrix& A);

/** 
 * Given a current error E = X - AD, dictionary D and coefficients A, update dictionary and error, 
 * so  that the total weight of the error E is reduced.
 * E = AD is an n x m matrix
 * D is a p x m matrix, where each row is an atom of dimension m
 * A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E
 */
idx_t update_dictionary_omp(binary_matrix& E, binary_matrix& D, const binary_matrix& A);

/**
 * dictionary update is done using PROXIMUS step
 */
idx_t update_dictionary_proximus_omp(binary_matrix& E, binary_matrix& D, binary_matrix& A);

/** 
 * Given a current error E, dictionary D and coefficients A, update coefficients and error, 
 * so that the total weight of the error E is reduced. The algorithm updates one row of A at a time.
 *
 * E = AD is an n x m matrix
 * D is a p x m matrix, where each row is an atom of dimension m
 * A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E 
*/
idx_t update_coefficients_omp(binary_matrix& E, const binary_matrix& D, binary_matrix& A);
#endif
