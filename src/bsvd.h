#ifndef BSVD_H
#define BSVD_H

#include "binmat.h"

void initialize_model_neighbor(const binary_matrix& E, 
		      binary_matrix& D, 
		      binary_matrix& A);

void initialize_model_random_centroids(const binary_matrix& E, 
		      binary_matrix& D, 
		      binary_matrix& A);

void initialize_model_random_centroids_xor(const binary_matrix& E, 
					   binary_matrix& D, 
					   binary_matrix& A);

void initialize_model_graph_grow(const binary_matrix& E, 
		      binary_matrix& D, 
		      binary_matrix& A);

void initialize_model_partition(const binary_matrix& E, 
		      binary_matrix& D, 
		      binary_matrix& A);

//#define initialize_model initialize_model_partition //  GOOD
//#define initialize_model initialize_model_neighbor  // BEST SO FAR
//#define initialize_model initialize_model_graph_grow VERY SLOW AND DOES NOT WORK WELL, DONT KNOW WHY
//#define initialize_model initialize_model_random_centroids WORST

/** 
 * Given a current error E = X - AD, dictionary D and coefficients A, update dictionary and error, 
 * so  that the total weight of the error E is reduced.
 * E = AD is an n x m matrix
 * D is a p x m matrix, where each row is an atom of dimension m
 * A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E
 */
idx_t update_dictionary_steepest(binary_matrix& E, binary_matrix& D, binary_matrix& A);
idx_t update_dictionary_steepest_omp(binary_matrix& E, binary_matrix& D, binary_matrix& A);

/**
 * dictionary update is done using PROXIMUS step
 */
idx_t update_dictionary_proximus(binary_matrix& E, binary_matrix& D, binary_matrix& A);
idx_t update_dictionary_proximus_omp(binary_matrix& E, binary_matrix& D, binary_matrix& A);

/** 
 * Given a current error E, dictionary D and coefficients A, update coefficients and error, 
 * so that the total weight of the error E is reduced. The algorithm updates one row of A at a time.
 *
 * E = AD is an n x m matrix
 * D is a p x m matrix, where each row is an atom of dimension m
 * A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E 
*/
idx_t update_coefficients_basic(binary_matrix& E, const binary_matrix& D, binary_matrix& A);

idx_t update_coefficients_omp(binary_matrix& E, const binary_matrix& D, binary_matrix& A);

/** 
 * Given a current error E, dictionary D and coefficients A, update coefficients and error, 
 * so that the total weight of the error E is reduced. The algorithm updates one row of A at a time.
 *
 * E = AD is an n x m matrix
 * D is a p x m matrix, where each row is an atom of dimension m
 * A is a n x p matrix, where each row contains the coefficients for representing the corresponding row of X=AD+E 
*/
idx_t update_coefficients_fast(binary_matrix& E, const binary_matrix& D, binary_matrix& A);

idx_t learn_model_traditional(binary_matrix& X,
			      binary_matrix& E, 
			      binary_matrix& D, 
			      binary_matrix& A);

idx_t learn_model_alter1(binary_matrix& X,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A);

idx_t learn_model_alter2(binary_matrix& X,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A);

idx_t learn_model_alter3(binary_matrix& X,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A);

idx_t learn_model_mdl_forward_selection(binary_matrix& X,
					binary_matrix& E, 
					binary_matrix& D, 
					binary_matrix& A);

idx_t learn_model_mdl_backward_selection(binary_matrix& X,
					 binary_matrix& E, 
					 binary_matrix& D, 
					 binary_matrix& A);

idx_t learn_model_mdl_full_search(binary_matrix& X,
				  binary_matrix& E, 
				  binary_matrix& D, 
				  binary_matrix& A);

typedef void (*mi_algorithm_t)(const binary_matrix& E, 
				binary_matrix& D, 
				binary_matrix& A);

typedef idx_t (*cu_algorithm_t)(binary_matrix& E, 
				const binary_matrix& D, 
				binary_matrix& A);

typedef idx_t (*du_algorithm_t)(binary_matrix& E, 
				binary_matrix& D, 
				binary_matrix& A);

typedef idx_t (*ml_algorithm_t)(binary_matrix& X,
				binary_matrix& E, 
				binary_matrix& D, 
				binary_matrix& A);

extern mi_algorithm_t initialize_model;
extern cu_algorithm_t update_coefficients;
extern du_algorithm_t update_dictionary;
extern ml_algorithm_t learn_model;
extern ml_algorithm_t learn_model_inner;

extern const char* mi_algorithm_names[];
extern const char* cu_algorithm_names[];
extern const char* du_algorithm_names[];
extern const char* lm_algorithm_names[];

extern long random_seed;

void learn_model_setup(int mi_algo, int cu_algo, int du_algo, int lm_algo, int lmi_algo);

#endif
