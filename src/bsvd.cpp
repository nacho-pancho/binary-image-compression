#include <algorithm>
#include "gsl/gsl_randist.h"
#include <omp.h>
#include "util.h"
#include <iomanip>
#include "bsvd.h"

static gsl_rng* get_rng() {
  static gsl_rng* rng = 0;
  if (rng == 0) { 
    rng = gsl_rng_alloc (gsl_rng_rand48);
    gsl_rng_set (rng, random_seed);                  // set random_seed
  } 
  return rng;
}

mi_algorithm_t initialize_model = initialize_model_neighbor;
cu_algorithm_t update_coefficients = update_coefficients_omp;
du_algorithm_t update_dictionary = update_dictionary;
ml_algorithm_t learn_model = learn_model_traditional;
ml_algorithm_t learn_model_inner = learn_model_traditional;

long random_seed = 34503498;

mi_algorithm_t mi_algorithm_catalog[] = {initialize_model_neighbor,
					 initialize_model_partition,
					 initialize_model_random_centroids,
					 initialize_model_random_centroids_xor,
					 initialize_model_graph_grow,
					 0};

const char* mi_algorithm_names[] = {"Neighbor initialization",
				    "Partition initialization",
				    "Random centroids initialization",
				    "Random centroids (in mod-2 algebra) initialization",
				    "Graph growing initialization",0
};

cu_algorithm_t cu_algorithm_catalog[] = {update_coefficients_omp,
					 update_coefficients,
					 update_coefficients_fast, // broken
					 0};

const char* cu_algorithm_names[] = {"OpenMP basic coefficients update",
				    "Basic coefficients update",
				    "Fast coefficients update (broken!)",0
};

du_algorithm_t du_algorithm_catalog[] = {update_dictionary_steepest,
					 update_dictionary_proximus,
					 update_dictionary_steepest_omp,
					 update_dictionary_proximus_omp,
					 0};

const char* du_algorithm_names[] = {"Steepest descent (a la MOD)  dictionary update",
				    "Proximus-like dictionary update",
				    "Steepest descent (a la MOD)  dictionary update (OMP)",
				    "Proximus-like dictionary update (OMP)",0
};

ml_algorithm_t learn_model_algorithm_catalog[] = {learn_model_traditional,
						  learn_model_alter1,
						  learn_model_alter2,
						  learn_model_alter3,
						  learn_model_mdl_forward_selection,
						  learn_model_mdl_backward_selection,
						  learn_model_mdl_full_search,
						  0};
const char* lm_algorithm_names[] = {"Model learning by traditional alternate descent",
				    "Role-switching learning 1: at each iteration, the role of A and D are switched",
				    "Role-switched learning 2: after convergence, the role of A and D are switched and traditional model is applied again",
				    "Role switched learning 3: like RS1 but only update_dictionary is applied (for use with Proximus",
				    "MDL/forward selection",
				    "MDO/backward selection",
				    "MDL/full search"
};


void learn_model_setup(int mi_algo, int cu_algo, int du_algo, int lm_algo, int lmi_algo) {
  if (mi_algo > 4) { std::cerr << "Invalid model initialization algorithm (0-" << 4 << ')' << std::endl; exit(-1); }
  if (cu_algo > 2) { std::cerr << "Invalid coefficients update algorithm (0-" << 2 << ')' << std::endl; exit(-1); }
  if (du_algo > 3) { std::cerr << "Invalid dictionary update algorithm (0-" << 3 << ')' << std::endl; exit(-1); }
  if (lm_algo > 6) { std::cerr << "Invalid model learning algorithm (0-" << 6 << ')' << std::endl; exit(-1); }
  if (lmi_algo > 3) { std::cerr << "Invalid inner model learning algorithm (0-" << 3 << ')' << std::endl; exit(-1); }

  initialize_model = mi_algorithm_catalog[mi_algo];
  std::cout << "Using " << mi_algorithm_names[mi_algo] << std::endl;
  update_coefficients = cu_algorithm_catalog[cu_algo];
  std::cout << "Using " << cu_algorithm_names[cu_algo] << std::endl;
  update_dictionary = du_algorithm_catalog[du_algo];
  std::cout << "Using " << du_algorithm_names[du_algo] << std::endl;
  learn_model = learn_model_algorithm_catalog[lm_algo];
  std::cout << "Using " << lm_algorithm_names[lm_algo] << " for outer learning loop." << std::endl;
  learn_model_inner = learn_model_algorithm_catalog[lmi_algo];
  std::cout << "Using " << lm_algorithm_names[lmi_algo] << " for inner learning." << std::endl;
}


void initialize_model_random_centroids_xor(const binary_matrix& E, 
					   binary_matrix& D, 
					   binary_matrix& A) {
  //
  // Initialize using random clusters
  //
  //
  // Initialize dictionary
  //
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix Ei(1,m);
  gsl_rng *rng = get_rng();  // random number generator
  A.clear();
  D.clear();
  for (idx_t i = 0; i < n; ++i) {
    idx_t k = (idx_t) gsl_rng_uniform_int(rng,p);
    A.set(i,k);
    D.copy_row_to(k,Dk);
    E.copy_row_to(i,Ei);
    add(Dk,Ei,Dk);
    D.set_row(k,Dk);
  }
  Ei.destroy();
  Dk.destroy();
}

void initialize_model_random_centroids(const binary_matrix& E, 
				       binary_matrix& D, 
				       binary_matrix& A) {
  //
  // Initialize using random clusters
  //
  //
  // Initialize dictionary
  //
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix Ei(1,m);
  gsl_rng *rng = get_rng();  // random number generator
  A.clear();
  D.clear();
  idx_t s[p][m];
  idx_t u[p];
  std::fill(&s[0][0],&s[0][0]+m*p,0);
  std::fill(u,u+p,0);
  idx_t S = 0;
  for (idx_t i = 0; i < n; ++i) {
    idx_t k = (idx_t) gsl_rng_uniform_int(rng,p);
    A.set(i,k);
    u[k]++;
    E.copy_row_to(i,Ei);
    for (idx_t j = 0; j < m; j++)
      if (Ei.get(0,j)) { S++; s[k][j]++; }
  }
  S /= (m*n);
  for (idx_t k = 0; k < p; k++) {
    for (idx_t j = 0; j < m; j++) {
      D.set(k,j, 2*s[k][j] >= u[k]);
    }
  }
  Ei.destroy();
  Dk.destroy();
}

/**
 * Works with up to m atoms in principle,
 * Each atom Dk is initialized to the centroid of all rows for which E
 * has ones on its k-th column
 */
void initialize_model_partition(const binary_matrix& E, 
				binary_matrix& D, 
				binary_matrix& A) {
  //
  // Initialize using random clusters
  //
  //
  // Initialize dictionary
  //
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix Ei(1,m);
  A.clear();
  D.clear();
  idx_t s[m];
  aux_t ranking[m];
  for (idx_t k = 0; k < m ; k++) {
    ranking[k].first = E.col_weight(k);
    ranking[k].second = k;
  }
  counting_sort(ranking,m);

  idx_t u;
  for (idx_t k = 0; (k < p) && (k < m); k++) {
    u = 0;
    std::fill(s,s+m,0);
    idx_t pivot = ranking[m-k-1].second;
    //    std::cout << "k=" << k << " pivot=" << pivot << " score=" << ranking[m-k-1].first << std::endl;
    for (idx_t i = 0; i < n; ++i) {
      if (E.get(i,pivot)) {
	u++;
	E.copy_row_to(i,Ei);
	for (idx_t j = 0; j < m; ++j) {
	  if (Ei.get(0,j))
	    s[j]++;
	}
      } // Ei had a row       
    }
    for (idx_t j = 0; j < m; ++j) 
      D.set(k,j,s[j] >= u/2);
  }
  // if p > m, the other atoms are left uninitialized!
  Ei.destroy();
  Dk.destroy();
}

/**
 * For each atom Dk we choose a random sample from E, Ej
 * and compute Dk as the centroid of all samples Ei that are neighbors to Ej in some way
 * in Boolean algebra, and PROXIMUS, this is interpreted as all other samples Ej that have a non-empty
 * intersection with Ei, that is, such that Ej AND Ei is non-zero. 
 */
void initialize_model_neighbor(const binary_matrix& E, 
			       binary_matrix& D, 
			       binary_matrix& A) {
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix Ei(1,m);
  binary_matrix Ej(1,m);
  gsl_rng *rng = get_rng();  // random number generator
  A.clear();
  D.clear();
  for (idx_t k = 0; k < p; ) {
    // pick a random row
    idx_t i = gsl_rng_uniform_int(rng,n);
    E.copy_row_to(i,Ei);
    if (Ei.weight()==0) continue;
    idx_t s[m];
    idx_t u = 0;
    std::fill(s,s+m,0);
    for (idx_t j = 0; j < n; ++j) {
      E.copy_row_to(j,Ej);
      bool_and(Ej,Ei,Ej);
      if (Ej.weight() > 0) {
	u++;
	for (idx_t j = 0; j < m; ++j) {
	  if (Ej.get(0,j))
	    s[j]++;
	}
      }
    }
    if (u > 0) { // the chosen has neighbors
      for (idx_t j = 0; j < m; ++j) 
	D.set(k,j,s[j] >= u/2);
      k++;
    }
  }
  Ei.destroy();
  Ej.destroy();
  Dk.destroy();
}

idx_t accumulate_to(binary_matrix& v, idx_t* s) {
  const idx_t m = v.get_cols();
  idx_t S=0;
  for (idx_t i = 0; i < m; i++) {
    if (v.get(0,i)) { 
      s[i]++;
      S++;
    }
  }
  return S;
}
/**
 * For each atom Dk we define a 'subgraph' and initialize it to a random row Ei, 
 * Then grow the subgraph by adding those rows Ej which share support with any element in the subgraph
 * I guess a good idea is to add not all, but
 * Finally, set Dk as the centroid of all samples in its corresponding part.
 */
void initialize_model_graph_grow(const binary_matrix& E, 
				 binary_matrix& D, 
				 binary_matrix& A) {
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix Ei(1,m);
  gsl_rng *rng = get_rng();  // random number generator
  A.clear();
  D.clear();
  idx_t s[p][m]; // part representation
  idx_t u[p]; // elements in each part
  bool t[n]; //mark if a given sample was selected
  std::fill(&s[0][0],&s[0][0]+m*p,0);
  std::fill(u,u+p,0);
  std::fill(t,t+n,false);
  //
  // initialize p 'subgraphs'. graphs are represented by a counts vector, which adds up 
  // the vectors in the graph.
  //
  int left = n;
  for (idx_t k = 0; (left >= 0) && (k < p); ) {
    // pick a random row
    idx_t i;
    do { 
      i = gsl_rng_uniform_int(rng,n);
    } while (t[i]);
    E.copy_row_to(i,Ei);
    for (idx_t j = 0; j <m; j++) {
      if (E.get(i,j)) s[k][j]=1; 
    }
    t[i] = true;
    left--;
    u[k] = 1;
    k++;
  }
  //
  // for each part, add to it the previously unused pattern that is closest
  // to its support. The original algorithm means ANY row that shares ANY non-zero 
  // with the part. HEre we pivot by weight: we prefer those that share the most used
  // atoms within that part.
  //
  while (left > 0) {
    for (idx_t k = 0; k < p; k++) {
      idx_t maxscore = 0;
      idx_t maxi = 0;
      idx_t score = 0;
      //
      // choose best newcomer
      //
      if (left <= 0) 
	break;
      for (idx_t i = 0; left && (i < n); i++) {
	if (t[i]) continue;
	for (idx_t j = 0; j <m; j++) {
	  //	  if (E.get(i,j)) score += s[k][j];  // NOT BAD< NEED TO TEST
	  if (E.get(i,j)) score  = 1; // do not accumulate overlapped supports, just support
	}
	if (score > maxscore) {
	  maxscore = score;
	  maxi = i;
	}
      }
      //std::cout << "left=" << left << " maxscore=" << maxscore << " maxi=" << maxi << " k=" << k << std::endl;
      if (maxscore == 0) { // reset part!
	idx_t i;
	do {
	  i = gsl_rng_uniform_int(rng,n);
	} while (t[i]);
	E.copy_row_to(i,Ei);
	for (idx_t j = 0; j < m; j++) {
	  s[k][j] = Ei.get(0,j) ? 1: 0;
	}
	t[i] = true;
	u[k] = 1;      
	left--;
      } else {
	//
	// add to part
	//
	t[maxi] = true;
	for (idx_t j = 0; j <m; j++) {
	  if (E.get(maxi,j)) s[k][j]++;
	}
	left--;
	u[k]++;
      }
    }
  }
  for (idx_t k = 0; k < p; k++) {
    for (idx_t j = 0; j < m; ++j) 
      D.set(k,j,s[k][j] >= u[k]/2);
  }
  
  Ei.destroy();
  Dk.destroy();
}
  
void initialize_model_random(const binary_matrix& E, 
			     binary_matrix& D, 
			     binary_matrix& A) {
  gsl_rng *rng = get_rng();  // random number generator
  const idx_t K = D.get_rows();
  const idx_t M = D.get_cols();
  for (idx_t k = 0; k < K; k++) {
    for (idx_t j = 0; j < M; j++) {
      D.set(k,j,gsl_ran_bernoulli (rng,0.5));
    }
  }
  A.clear();
}

idx_t update_coefficients_basic(binary_matrix& E, const binary_matrix& D,  binary_matrix& A) 
{
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  std::cout << "cu/basic" << std::endl;
  //
  // SPARSE CODING STEP
  // Go over each row and encode it using rows from D until the residual weight
  // is no longer diminished.
  //
  idx_t changed = 0;
  binary_matrix Ei(1,m);
  binary_matrix Ai(1,p);
  binary_matrix Dk(1,m);
  for (idx_t i = 0; i < n; i++) {
    E.copy_row_to(i,Ei);
    A.copy_row_to(i,Ai);
    bool improved = true;
    bool ichanged = false;
    idx_t iter = 0;
    while (improved) {
      idx_t w = Ei.weight();
      D.copy_row_to(0,Dk);
      idx_t bestk = 0, bestd = dist(Ei,Dk);
      bool_and(Ei,Dk,Dk);
      //     idx_t r = Dk.weight();
      //      std::cout << "\titer=" << iter << " k=" << 0 << " r=" << r << " d=" << bestd << std::endl;
      for (idx_t k = 1; k < p; k++) {
	D.copy_row_to(k,Dk);
	const idx_t dk = dist(Ei,Dk);
	bool_and(Ei,Dk,Dk);
	//	idx_t r = Dk.weight();
	//	std::cout << "\titer" << iter << " k=" << k << " r=" << r << " d=" << dk << std::endl;
	if (dk < bestd) {
	  bestd = dk;
	  bestk = k;
	} 
      }
      //      std::cout << "i=" << i << " w=" << w << " bestk=" << bestk << " bestd=" << bestd << std::endl;
      if (bestd < w) {
	D.copy_row_to(bestk,Dk);
	Ai.flip(0,bestk);
	add(Ei,Dk,Ei);
	w = bestd;
	ichanged = true;
      } else {
	improved = false;
      }  
      iter++;
    } // while there is any improvement
    if (ichanged) {
      changed++;
      E.set_row(i,Ei);
      A.set_row(i,Ai);
    }
  }
  Ei.destroy();
  Ai.destroy();
  Dk.destroy();
  return changed;
}


idx_t update_dictionary_steepest(binary_matrix& E, binary_matrix& D,  binary_matrix& A)
{
  //
  // DICTIONARY UPDATE STEP
  // For each atom,
  // take the samples that use it
  // if for the j-th dimension the weight is larger than m/2, turn j-th
  // atom element on, else turn it off
  //
  //  std::cout << "du/stee" << std::endl;
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix newDk(1,m);
  binary_matrix Ei(1,m);

  idx_t changed = 0;
  for (idx_t k = 0; k < p; k++) {
    D.copy_row_to(k,Dk);
    idx_t weights[(int)m];
    idx_t atom_usage = 0;
    std::fill(weights,weights+m,0);
    for (idx_t i = 0; i < n; i++) {
      idx_t aik = A.get(i,k);
      if (aik) { // if atom k was used in i-th samples
	atom_usage++;
	E.copy_row_to(i,Ei);
	add(Ei,Dk,Ei); // add-back old atom, we do not take it into account
	for (idx_t j = 0; j < m ; j++) {
	  if (Ei.get(0,j)) {
	    weights[j]++;
	  }
	} // add i-th sample to weights vector
      } // if atom was used
    } // for each data sample
    if (!atom_usage) 
      continue;
    // now update atom based on weights
    const idx_t u = atom_usage/2;
    Dk.copy_to(newDk);
    for (idx_t j = 0; j < m ; j++) {
      newDk.set(0,j, weights[j] > u);
    }
    if (dist(newDk,Dk) > 0) { 
      // there was a change in the atom
      changed++;
      D.set_row(k,newDk);
      // update residual E
      for (idx_t i = 0; i < n; i++) {
	idx_t aik = A.get(i,k);
	if (!aik) 
	  continue;
	E.copy_row_to(i,Ei);
	add(Ei,Dk,Ei); // add-back old atom
	add(Ei,newDk,Ei); // substract (same as add here) new atom
	E.set_row(i,Ei);
      }
    }
  } // for each atom
  newDk.destroy();
  Dk.destroy();
  Ei.destroy();
  return changed;
} // end


idx_t update_dictionary_proximus(binary_matrix& E, binary_matrix& D, binary_matrix& A)
{
  //
  // DICTIONARY UPDATE STEP
  // For each atom dk,
  // take the subset of samples that use it
  // remove dk from them
  // apply the proximus rank-one approximation algorithm to dk and a^k
  // 
  //  std::cout << "du/prox" << std::endl;
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Ei(1,m); // rows of E
  binary_matrix Ej(1,n); // column of E
  binary_matrix Dk(1,m); // column of D
  binary_matrix Ak(1,n); // row of coefs
  binary_matrix newDk(1,m); // column of D
  binary_matrix newAk(1,n); // row of coefs


  idx_t changed = 0;
  //
  // FOR EACH ATOM
  //
  for (idx_t k = 0; k < p; k++) {
    //
    // PROXIMUS-LIKE ITERATION, UNTIL NO CHANGE OCCURS
    // 
    bool kchanged = false;
    bool converged;
    //    std::cout << "PROXIMUS: k=" << k << std::endl;
    //
    // PROXIMUS INITIALIZATION of Ak: maximize correlation, not distance
    //
    // ----------------
#if 0
    A.copy_col_to(k,Ak);
    idx_t da = 0;
    //
    // i. compute s = <E,D_k^t>
    //
    std::pair<idx_t,idx_t> s[(int)n];
    for (idx_t i = 0 ; i < n; i++) {
      s[i].first = 0; s[i].second = i;
    }
    for (idx_t j = 0; j < m; j++) {
      idx_t dkj = D.get(k,j);
      if (dkj) { // if atom k was used in i-th samples
	E.copy_col_to(j,Ej);
	add(Ej,Ak,Ej); // add-back old coefs, we do not take it into account
	for (idx_t i = 0; i < n ; i++) {
	  if (Ej.get(0,i)) {
	    s[i].first++;
	  }
	} // add i-th sample to weights vector
      } // if atom was used
    } // for each data sample
      //
      // ii. find A_k that maximizes |A_k^t s| / |A_k|
      // this is simply setting the p largest elements of E D_k^t, where p is such that 
      // z_(p+1) < \sum_{k=1}^{p} z_k / p, where z is the sorted version of s
      // THIS SHOULD BE REPLACED WITH COUNTING SORT, WHICH IS O(n)
      // O(nlog n) kills the algorithm
    std::sort(s,s+n);    
    std::cout << s[n-1].first << std::endl;
    //counting_sort(s,n);
    newAk.clear();
    for (int i = (n-1), sp = 0; (i>=0) && (s[i].first > sp); i--) {
      newAk.set(0,s[i].second);
      sp += s[i].first;
    }
    da = dist(newAk,Ak);
    if (da > 0) { 
      // there was a change 
      A.set_col(k,newAk);
      // update residual E
      for (idx_t j = 0; j < m; j++) {
	idx_t dkj = D.get(k,j);
	if (!dkj) 
	  continue;
	E.copy_col_to(j,Ej);
	add(Ej,Ak,Ej); // add-back old atom
	add(Ej,newAk,Ej); // substract (same as add here) new atom
	E.set_col(j,Ej);
      }
    } else {
      continue;
    }
    
#endif
    // ----------------

    do { // until converged is true
      converged = true;
      //
      // Update atom Dk
      //
      idx_t u = 0;
      idx_t Dw[(int)m];
      std::fill(Dw,Dw+m,0);
      D.copy_row_to(k,Dk);

      for (idx_t i = 0; i < n; i++) {
	idx_t aik = A.get(i,k);
	if (aik) { // if atom k was used in i-th samples
	  u++;
	  E.copy_row_to(i,Ei);
	  add(Ei,Dk,Ei); // add-back old atom, we do not take it into account
	  for (idx_t j = 0; j < m ; j++) {
	    if (Ei.get(0,j)) {
	      Dw[j]++;
	    }
	  } // add i-th sample to weights vector
	} // if atom was used
      } // for each data sample
      idx_t dd = 0;
      if (u) { 
	// now update atom based on weights
	u /= 2; 
	Dk.copy_to(newDk);
	for (idx_t j = 0; j < m ; j++) {
	  newDk.set(0,j, Dw[j] > u);
	}
	dd = dist(newDk,Dk);
	if (dd > 0) { 
	  // there was a change in the atom
	  D.set_row(k,newDk);
	  converged = false;
	  kchanged = true;
	  // update residual matrix E
	  for (idx_t i = 0; i < n; i++) {
	    idx_t aik = A.get(i,k);
	    if (!aik) 
	      continue;
	    E.copy_row_to(i,Ei);
	    add(Ei,Dk,Ei); // add-back old atom
	    add(Ei,newDk,Ei); // substract (same as add here) new atom
	    E.set_row(i,Ei);
	  }
	}
      }
      //      std::cout << "\tu=" << u << "\tdd=" << dd << "\t|E|=" << E.weight() << std::endl;      
      //
      // UPDATE COEFs, Ak
      //
      u = 0;
      A.copy_col_to(k,Ak);
      idx_t Aw[(int)n];
      std::fill(Aw,Aw+n,0);
      idx_t da = 0;
      for (idx_t j = 0; j < m; j++) {
	idx_t dkj = D.get(k,j);
	if (dkj) { // if atom k was used in i-th samples
	  u++;
	  E.copy_col_to(j,Ej);
	  add(Ej,Ak,Ej); // add-back old coefs, we do not take it into account
	  for (idx_t i = 0; i < n ; i++) {
	    if (Ej.get(0,i)) {
	      Aw[i]++;
	    }
	  } // add i-th sample to weights vector
	} // if atom was used
      } // for each data sample
      if (u) {
	// now update atom based on weights
	u /= 2;
	Ak.copy_to(newAk);
	for (idx_t i = 0; i < n ; i++) {
	  newAk.set(0,i, Aw[i] > u);
	}
	da = dist(newAk,Ak);
	if (da > 0) { 
	  // there was a change 
	  A.set_col(k,newAk);
	  converged = false;
	  //	  kchanged = true;
	  // update residual E
	  for (idx_t j = 0; j < m; j++) {
	    idx_t dkj = D.get(k,j);
	    if (!dkj) 
	      continue;
	    E.copy_col_to(j,Ej);
	    add(Ej,Ak,Ej); // add-back old atom
	    add(Ej,newAk,Ej); // substract (same as add here) new atom
	    E.set_col(j,Ej);
	  }
	}
	//	std::cout << "\tu=" << u << "\tda=" << da << "\t|E|=" << E.weight() << std::endl;      
      }
    } while (!converged);
    //
    // the pair Dk, Ak has changed during this run
    //
    if (kchanged)
      changed++;
  } // for each atom
    // finished
  Ei.destroy();
  Ej.destroy();
  newDk.destroy();
  newAk.destroy();
  Dk.destroy();
  Ak.destroy();
  return changed;
} // end


idx_t update_dictionary_steepest_omp(binary_matrix& E, binary_matrix& D,  binary_matrix& A)
{
  //
  // DICTIONARY UPDATE STEP
  // For each atom,
  // take the samples that use it
  // if for the j-th dimension the weight is larger than m/2, turn j-th
  // atom element on, else turn it off
  //
  //  std::cout << "du/stee" << std::endl;
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Dk(1,m);
  binary_matrix newDk(1,m);

  omp_set_num_threads(omp_get_max_threads());
  idx_t NT;
#pragma omp parallel
  {  
 NT = omp_get_num_threads(); 
  }
  //  std::cout << "THREADS=" << NT << std::endl;
  binary_matrix Ei[NT];
  for (idx_t TT = 0; TT < NT; TT++) {
    Ei[TT].allocate(1,m);
  }
  const binary_matrix cA = A;
  //  binary_matrix Ei(1,m);

  idx_t changed = 0;
  for (idx_t k = 0; k < p; k++) {
    const idx_t T = omp_get_thread_num();    
    D.copy_row_to(k,Dk);
    idx_t weights[(int)m];
    idx_t atom_usage = 0;
    std::fill(weights,weights+m,0);
#pragma omp parallel for 
    for (idx_t i = 0; i < n; i++) {
      if (cA.get(i,k)) { // if atom k was used in i-th samples
	atom_usage++;
	E.copy_row_to(i,Ei[T]);
	add(Ei[T],Dk,Ei[T]); // add-back old atom, we do not take it into account
	for (idx_t j = 0; j < m ; j++) {
	  if (Ei[T].get(0,j)) {
	    weights[j]++;
	  }
	} // add i-th sample to weights vector
      } // if atom was used
    } // for each data sample

    if (!atom_usage) 
      continue;
    // now update atom based on weights
    const idx_t u = atom_usage/2;
    Dk.copy_to(newDk);
    for (idx_t j = 0; j < m ; j++) {
      newDk.set(0,j, weights[j] > u);
    }
    if (dist(newDk,Dk) > 0) { 
      // there was a change in the atom
      changed++;
      D.set_row(k,newDk);
      // update residual E
#pragma omp parallel for 
      for (idx_t i = 0; i < n; i++) {
	if (!cA.get(i,k)) 
	  continue;
	E.copy_row_to(i,Ei[T]);
	add(Ei[T],Dk,Ei[T]); // add-back old atom
	add(Ei[T],newDk,Ei[T]); // substract (same as add here) new atom
	E.set_row(i,Ei[T]);
      }
    }
  } // for each atom
  newDk.destroy();
  Dk.destroy();
  for (idx_t T = 0; T < NT; T++) {
    Ei[T].destroy();
  }
  return changed;
} // end


idx_t update_dictionary_proximus_omp(binary_matrix& E, binary_matrix& D, binary_matrix& A)
{
  //
  // DICTIONARY UPDATE STEP
  // For each atom dk,
  // take the subset of samples that use it
  // remove dk from them
  // apply the proximus rank-one approximation algorithm to dk and a^k
  // 
  //  std::cout << "du/prox" << std::endl;
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  binary_matrix Ei(1,m); // rows of E
  binary_matrix Ej(1,n); // column of E
  binary_matrix Dk(1,m); // column of D
  binary_matrix Ak(1,n); // row of coefs
  binary_matrix newDk(1,m); // column of D
  binary_matrix newAk(1,n); // row of coefs


  idx_t changed = 0;
  //
  // FOR EACH ATOM
  //
  for (idx_t k = 0; k < p; k++) {
    //
    // PROXIMUS-LIKE ITERATION, UNTIL NO CHANGE OCCURS
    // 
    bool kchanged = false;
    bool converged;
    //    std::cout << "PROXIMUS: k=" << k << std::endl;
    //
    // PROXIMUS INITIALIZATION of Ak: maximize correlation, not distance
    //
    // ----------------
#if 0
    A.copy_col_to(k,Ak);
    idx_t da = 0;
    //
    // i. compute s = <E,D_k^t>
    //
    std::pair<idx_t,idx_t> s[(int)n];
    for (idx_t i = 0 ; i < n; i++) {
      s[i].first = 0; s[i].second = i;
    }
    for (idx_t j = 0; j < m; j++) {
      idx_t dkj = D.get(k,j);
      if (dkj) { // if atom k was used in i-th samples
	E.copy_col_to(j,Ej);
	add(Ej,Ak,Ej); // add-back old coefs, we do not take it into account
	for (idx_t i = 0; i < n ; i++) {
	  if (Ej.get(0,i)) {
	    s[i].first++;
	  }
	} // add i-th sample to weights vector
      } // if atom was used
    } // for each data sample
      //
      // ii. find A_k that maximizes |A_k^t s| / |A_k|
      // this is simply setting the p largest elements of E D_k^t, where p is such that 
      // z_(p+1) < \sum_{k=1}^{p} z_k / p, where z is the sorted version of s
      // THIS SHOULD BE REPLACED WITH COUNTING SORT, WHICH IS O(n)
      // O(nlog n) kills the algorithm
    std::sort(s,s+n);    
    std::cout << s[n-1].first << std::endl;
    //counting_sort(s,n);
    newAk.clear();
    for (int i = (n-1), sp = 0; (i>=0) && (s[i].first > sp); i--) {
      newAk.set(0,s[i].second);
      sp += s[i].first;
    }
    da = dist(newAk,Ak);
    if (da > 0) { 
      // there was a change 
      A.set_col(k,newAk);
      // update residual E
      for (idx_t j = 0; j < m; j++) {
	idx_t dkj = D.get(k,j);
	if (!dkj) 
	  continue;
	E.copy_col_to(j,Ej);
	add(Ej,Ak,Ej); // add-back old atom
	add(Ej,newAk,Ej); // substract (same as add here) new atom
	E.set_col(j,Ej);
      }
    } else {
      continue;
    }
    
#endif
    // ----------------

    do { // until converged is true
      converged = true;
      //
      // Update atom Dk
      //
      idx_t u = 0;
      idx_t Dw[(int)m];
      std::fill(Dw,Dw+m,0);
      D.copy_row_to(k,Dk);

      for (idx_t i = 0; i < n; i++) {
	idx_t aik = A.get(i,k);
	if (aik) { // if atom k was used in i-th samples
	  u++;
	  E.copy_row_to(i,Ei);
	  add(Ei,Dk,Ei); // add-back old atom, we do not take it into account
	  for (idx_t j = 0; j < m ; j++) {
	    if (Ei.get(0,j)) {
	      Dw[j]++;
	    }
	  } // add i-th sample to weights vector
	} // if atom was used
      } // for each data sample
      idx_t dd = 0;
      if (u) { 
	// now update atom based on weights
	u /= 2; 
	Dk.copy_to(newDk);
	for (idx_t j = 0; j < m ; j++) {
	  newDk.set(0,j, Dw[j] > u);
	}
	dd = dist(newDk,Dk);
	if (dd > 0) { 
	  // there was a change in the atom
	  D.set_row(k,newDk);
	  converged = false;
	  kchanged = true;
	  // update residual matrix E
	  for (idx_t i = 0; i < n; i++) {
	    idx_t aik = A.get(i,k);
	    if (!aik) 
	      continue;
	    E.copy_row_to(i,Ei);
	    add(Ei,Dk,Ei); // add-back old atom
	    add(Ei,newDk,Ei); // substract (same as add here) new atom
	    E.set_row(i,Ei);
	  }
	}
      }
      //      std::cout << "\tu=" << u << "\tdd=" << dd << "\t|E|=" << E.weight() << std::endl;      
      //
      // UPDATE COEFs, Ak
      //
      u = 0;
      A.copy_col_to(k,Ak);
      idx_t Aw[(int)n];
      std::fill(Aw,Aw+n,0);
      idx_t da = 0;
      for (idx_t j = 0; j < m; j++) {
	idx_t dkj = D.get(k,j);
	if (dkj) { // if atom k was used in i-th samples
	  u++;
	  E.copy_col_to(j,Ej);
	  add(Ej,Ak,Ej); // add-back old coefs, we do not take it into account
	  for (idx_t i = 0; i < n ; i++) {
	    if (Ej.get(0,i)) {
	      Aw[i]++;
	    }
	  } // add i-th sample to weights vector
	} // if atom was used
      } // for each data sample
      if (u) {
	// now update atom based on weights
	u /= 2;
	Ak.copy_to(newAk);
	for (idx_t i = 0; i < n ; i++) {
	  newAk.set(0,i, Aw[i] > u);
	}
	da = dist(newAk,Ak);
	if (da > 0) { 
	  // there was a change 
	  A.set_col(k,newAk);
	  converged = false;
	  //	  kchanged = true;
	  // update residual E
	  for (idx_t j = 0; j < m; j++) {
	    idx_t dkj = D.get(k,j);
	    if (!dkj) 
	      continue;
	    E.copy_col_to(j,Ej);
	    add(Ej,Ak,Ej); // add-back old atom
	    add(Ej,newAk,Ej); // substract (same as add here) new atom
	    E.set_col(j,Ej);
	  }
	}
	//	std::cout << "\tu=" << u << "\tda=" << da << "\t|E|=" << E.weight() << std::endl;      
      }
    } while (!converged);
    //
    // the pair Dk, Ak has changed during this run
    //
    if (kchanged)
      changed++;
  } // for each atom
    // finished
  Ei.destroy();
  Ej.destroy();
  newDk.destroy();
  newAk.destroy();
  Dk.destroy();
  Ak.destroy();
  return changed;
} // end

idx_t update_coefficients_omp(binary_matrix& E, const binary_matrix& D,  binary_matrix& A) 
{
  //  std::cout << "cu/omp" << std::endl;
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  //
  // SPARSE CODING STEP
  // Go over each row and encode it using rows from D until the residual weight
  // is no longer diminished.
  //
  idx_t changed = 0;
  omp_set_num_threads(omp_get_max_threads());
  idx_t NT;
#pragma omp parallel
  {  
 NT = omp_get_num_threads(); 
  }
  //  std::cout << "THREADS=" << NT << std::endl;
  binary_matrix Ei[NT];
  binary_matrix Ai[NT];
  binary_matrix Dk[NT];
  for (idx_t TT = 0; TT < NT; TT++) {
    Ei[TT].allocate(1,m);
    Ai[TT].allocate(1,p);
    Dk[TT].allocate(1,m);
  }

#pragma omp parallel for schedule(dynamic)
  for (idx_t i = 0; i < n; i++) {
    const idx_t T = omp_get_thread_num();
    E.copy_row_to(i,Ei[T]);
    A.copy_row_to(i,Ai[T]);
    bool improved = true;
    bool ichanged = false;
    idx_t iter = 0;
    while (improved) {
      idx_t w = Ei[T].weight();
      D.copy_row_to(0,Dk[T]);
      idx_t bestk = 0, bestd = dist(Ei[T],Dk[T]);
      bool_and(Ei[T],Dk[T],Dk[T]);
      //     idx_t r = Dk[T].weight();
      //      std::cout << "\titer=" << iter << " k=" << 0 << " r=" << r << " d=" << bestd << std::endl;
      for (idx_t k = 1; k < p; k++) {
	D.copy_row_to(k,Dk[T]);
	const idx_t dk = dist(Ei[T],Dk[T]);
	bool_and(Ei[T],Dk[T],Dk[T]);
	//	idx_t r = Dk[T].weight();
	//	std::cout << "\titer" << iter << " k=" << k << " r=" << r << " d=" << dk << std::endl;
	if (dk < bestd) {
          bestd = dk;
          bestk = k;
	} 
      }
      //      std::cout << "i=" << i << " w=" << w << " bestk=" << bestk << " bestd=" << bestd << std::endl;
      if (bestd < w) {
	D.copy_row_to(bestk,Dk[T]);
	Ai[T].flip(0,bestk);
	add(Ei[T],Dk[T],Ei[T]);
	w = bestd;
	ichanged = true;
      } else {
	improved = false;
      }  
      iter++;
    } // while there is any improvement
    if (ichanged) {
      changed++;
      E.set_row(i,Ei[T]);
      A.set_row(i,Ai[T]);
    }
  }
  for (idx_t T = 0; T < NT; T++) {
    Ei[T].destroy();
    Ai[T].destroy();
    Dk[T].destroy();
  }
  return changed;
}


  //
  // DOES NOT WORK WELL!
  // I need to find proper surrogates to the distance between E and D that
  // can actually be stored.
  // But speed is of no concern now.
  //
idx_t update_coefficients_fast(binary_matrix& E, const binary_matrix& D,  binary_matrix& A) 
{
  //  std::cout << "cu/fast" << std::endl;
  const idx_t m = E.get_cols();
  const idx_t n = E.get_rows();
  const idx_t p = D.get_rows();
  //
  // SPARSE CODING STEP
  // Go over each row and encode it using rows from D until the residual weight
  // is no longer diminished.
  //
  idx_t changed = 0;
  idx_t G[p][p];
  idx_t r[p];
  idx_t h[p];
  binary_matrix dk(1,m),eik(1,m);
  binary_matrix& d1 = dk; // aliases
  binary_matrix& d2 = eik;
  for (idx_t k1 = 0; k1 < p; ++k1) {
    D.copy_row_to(k1,d1);
    for (idx_t k2 = k1; k2 < p; ++k2) {
      D.copy_row_to(k2,d2);
      bool_and(d1,d2,d2);
      G[k1][k2] = G[k2][k1] = d2.weight();
    }
  }
#if 0
  for (idx_t k1 = 0; k1 < p; ++k1) {
    for (idx_t k2 = 0; k2 < p; ++k2) {
      std::cout << " k1=" << k1 << " k2=" << k2 << " G=" << G[k1][k2] << std::endl;
    }
  }
#endif
  binary_matrix Ei(1,m);
  binary_matrix Ai(1,p); 
  for (idx_t i = 0; i < n; i++) {
    E.copy_row_to(i,Ei);
    A.copy_row_to(i,Ai);
    idx_t w = Ei.weight();
    // EitD^(0)
    // h(0)
    bool ichanged = false;
    idx_t iter = 0;
    for (idx_t k = 0; k < p; ++k) {
      D.copy_row_to(k,dk);
      bool_and(dk,Ei,eik);
      r[k] = eik.weight();
      h[k] = G[k][k] - 2*r[k] + w;
      //std::cout << "\titer=" << iter << " k=" << k1 << " r=" << r[k1] << " d=" << h[k1] << std::endl;
    }
    // a(0) = 0
    bool improved = true;
    while (improved) {
      const idx_t bestk = std::min_element(h,h+p) - h;
      const idx_t bestd = h[bestk];
      //std::cout << "i=" << i << " w=" << w << " bestk=" << bestk << " bestd=" << bestd << std::endl;
      if (bestd < w) {
	idx_t olda = Ai.get(0,bestk);
	Ai.flip(0,bestk);
	D.copy_row_to(bestk,dk);
	add(Ei,dk,Ei);
	w = bestd;
	idx_t oldrbestk = r[bestk];
	iter++;
	if (olda) { // remove atom from solution
	  for (idx_t k = 0; k < p; ++k) {
	    //std::cout << "G[k][bestk]=" << G[k][bestk] << std::endl;
	    r[k] = r[k] + G[k][bestk];  
	    h[k] = h[k] + G[bestk][bestk] + 2*(((int)oldrbestk-G[k][bestk]));  
	    D.copy_row_to(k,dk); // DEBUG!
	    //std::cout << "\titer=" << iter << " k=" << k << " r=" << r[k] << " d=" << h[k] << " d=" << dist(Ei,dk) << std::endl;
	  }
	} else  { // add atom to solution
	  for (idx_t k = 0; k < p; ++k) {
	    //std::cout << "G[k][bestk]=" << G[k][bestk] << std::endl;
	    r[k] = r[k] - G[k][bestk];  
	    h[k] = h[k] + G[bestk][bestk] - 2*(((int)oldrbestk-G[k][bestk]));  
	    D.copy_row_to(k,dk); // DEBUG!
	    //std::cout << "\titer=" << iter << " k=" << k << " r=" << r[k] << " d=" << h[k] << " d=" << dist(Ei,dk) << std::endl;
	  }
	}
	ichanged = true;
      } else {
	improved = false;
      } 
    } // while there is any improvement
    if (ichanged) {
      changed++;
      E.set_row(i,Ei);
      A.set_row(i,Ai);
    }
  }
  dk.destroy();
  eik.destroy();
  Ei.destroy();
  Ai.destroy();
  return changed;
}

idx_t learn_model_traditional(binary_matrix& X,
			      binary_matrix& E, 
			      binary_matrix& D, 
			      binary_matrix& A) {
  mul(A,false,D,false,E);
  add(E,X,E);
  //  std::cout << "trad" << std::endl;
  idx_t changed = 1;
  idx_t iter = 0;
  //std::cout << "iter=" << std::setw(8) << iter 
//	    << "\t||E||=" << std::setw(8) << E.weight()
//	    << "\t||D||=" << std::setw(8) << D.weight()
//	    << "\t||A||=" << std::setw(8) << A.weight() << std::endl;
  while (changed > 0) {    
    iter++;
    idx_t changed_coefs = update_coefficients(E,D,A);
  //  std::cout << "iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
    changed = changed_coefs + update_dictionary(E,D,A);
//    std::cout << "iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;
  }
  return iter;
}


idx_t learn_model_alter1(binary_matrix& X,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A) {
  //  std::cout << "alter1" << std::endl;
  const idx_t N = E.get_rows();
  const idx_t M = E.get_cols();
  const idx_t K = D.get_rows();

  mul(A,false,D,false,E);
  add(E,X,E);
  binary_matrix Dt(M,K);
  binary_matrix At(K,N);
  binary_matrix Et(M,N);

  //
  // RUN BSVD
  //
  idx_t changed = 1;
  idx_t iter = 0;
//  std::cout << "iter=" << std::setw(8) << iter 
//	    << "\t||E||=" << std::setw(8) << E.weight()
//	    << "\t||D||=" << std::setw(8) << D.weight()
//	    << "\t||A||=" << std::setw(8) << A.weight() << std::endl;
  while (changed > 0) {    
    iter++;
    idx_t changed_coefs = update_coefficients(E,D,A);
//    std::cout << "DIRECT: iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
    changed = changed_coefs + update_dictionary(E,D,A);
  //  std::cout << "DIRECT: iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;

    A.transpose_to(At);
    D.transpose_to(Dt);
    E.transpose_to(Et);

    changed_coefs = update_coefficients(Et,At,Dt);
  //  std::cout << "TRANSP: iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << Et.weight()
//	      << "\t||D||=" << std::setw(8) << Dt.weight()
//	      << "\t||A||=" << std::setw(8) << At.weight()
//	      << "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;

    changed = update_dictionary(Et,At,Dt);
  //  std::cout << "TRANSP: iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << Et.weight()
//	      << "\t||D||=" << std::setw(8) << Dt.weight()
//	      << "\t||A||=" << std::setw(8) << At.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;
    
    At.transpose_to(A);
    Dt.transpose_to(D);
    Et.transpose_to(E);
  }
  At.destroy();
  Dt.destroy();
  Et.destroy();
  return iter;
}


idx_t learn_model_alter2(binary_matrix& X,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A) {
  //  std::cout << "alter2" << std::endl;
  const idx_t N = E.get_rows();
  const idx_t M = E.get_cols();
  const idx_t K = D.get_rows();

  mul(A,false,D,false,E);
  add(E,X,E);
  binary_matrix Dt(M,K);
  binary_matrix At(K,N);
  binary_matrix Et(M,N);
  idx_t changed = 1;
  idx_t iter = 0;
//  std::cout << "iter=" << std::setw(8) << iter 
//	    << "\t||E||=" << std::setw(8) << E.weight()
//	    << "\t||D||=" << std::setw(8) << D.weight()
//	    << "\t||A||=" << std::setw(8) << A.weight() << std::endl;

  idx_t outer_changed = 1;
  while (outer_changed > 0) {
    outer_changed = 0;
    while (changed > 0) {    
      iter++;
      idx_t changed_coefs = update_coefficients(E,D,A);
//      std::cout << "DIRECT: iter=" << std::setw(8) << iter 
//		<< "\t||E||=" << std::setw(8) << E.weight()
//		<< "\t||D||=" << std::setw(8) << D.weight()
//		<< "\t||A||=" << std::setw(8) << A.weight()
//		<< "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;
      changed = changed_coefs + update_dictionary(E,D,A);
//      std::cout << "DIRECT: iter=" << std::setw(8) << iter 
//		<< "\t||E||=" << std::setw(8) << E.weight()
//		<< "\t||D||=" << std::setw(8) << D.weight()
//		<< "\t||A||=" << std::setw(8) << A.weight()
//		<< "\tchanged atoms=" << std::setw(8) << changed << std::endl;
      outer_changed += changed;
    }
    A.transpose_to(At);
    D.transpose_to(Dt);
    E.transpose_to(Et);
    changed = 1;
    iter = 0;
    while (changed > 0) {    
      iter++;
      idx_t changed_coefs = update_coefficients(Et,At,Dt);
//      std::cout << "TRANSPOSED: iter=" << std::setw(8) << iter 
//		<< "\t||E||=" << std::setw(8) << Et.weight()
//		<< "\t||D||=" << std::setw(8) << Dt.weight()
//		<< "\t||A||=" << std::setw(8) << At.weight()
//		<< "\tchanged coefs=" << std::setw(8) << changed_coefs << std::endl;

      changed = changed_coefs + update_dictionary(Et,At,Dt);
//      std::cout << "TRANSPOSED: iter=" << std::setw(8) << iter 
//		<< "\t||E||=" << std::setw(8) << Et.weight()
//		<< "\t||D||=" << std::setw(8) << Dt.weight()
//		<< "\t||A||=" << std::setw(8) << At.weight()
//		<< "\tchanged atoms=" << std::setw(8) << changed << std::endl; 
      outer_changed += changed;   
    }
    At.transpose_to(A);
    Dt.transpose_to(D);
    Et.transpose_to(E);
  }
  At.destroy();
  Dt.destroy();
  Et.destroy();
  return iter;
}


idx_t learn_model_alter3(binary_matrix& X,
			 binary_matrix& E, 
			 binary_matrix& D, 
			 binary_matrix& A) {
  //  std::cout << "alter3" << std::endl;
  const idx_t N = E.get_rows();
  const idx_t M = E.get_cols();
  const idx_t K = D.get_rows();
  mul(A,false,D,false,E);
  add(E,X,E);

  binary_matrix Dt(M,K);
  binary_matrix At(K,N);
  binary_matrix Et(M,N);
  idx_t changed = K+1;
  idx_t iter = 0;
//  std::cout << "iter=" << std::setw(8) << iter 
//	    << "\t||E||=" << std::setw(8) << E.weight()
//	    << "\t||D||=" << std::setw(8) << D.weight()
//	    << "\t||A||=" << std::setw(8) << A.weight() << std::endl;
  while (changed > 0) {    
    iter++;
    A.transpose_to(At);
    D.transpose_to(Dt);
    E.transpose_to(Et);
    changed = update_dictionary(Et,At,Dt);
//    std::cout << "iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << Et.weight()
//	      << "\t||D||=" << std::setw(8) << Dt.weight()
//	      << "\t||A||=" << std::setw(8) << At.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;

    At.transpose_to(A);
    Dt.transpose_to(D);
    Et.transpose_to(E);
    changed = update_dictionary(E,D,A);
//    std::cout << "iter=" << std::setw(8) << iter 
//	      << "\t||E||=" << std::setw(8) << E.weight()
//	      << "\t||D||=" << std::setw(8) << D.weight()
//	      << "\t||A||=" << std::setw(8) << A.weight()
//	      << "\tchanged atoms=" << std::setw(8) << changed << std::endl;
  }
  At.destroy();
  Dt.destroy();
  Et.destroy();
return iter;
}

#include "coding.h"

idx_t model_codelength(const binary_matrix& E, 
		       const binary_matrix& D, 
		       const binary_matrix& A) {

  const idx_t M = E.get_cols();
  const idx_t N = E.get_rows();
  const idx_t K = D.get_rows();


  binary_matrix Dk(1,M);
  binary_matrix Ak(1,N);

  idx_t LE = universal_codelength(E.get_rows()*E.get_cols(),E.weight());
  idx_t LD = 0, LA = 0;
  for (idx_t k = 0; k < K; k++) {
    D.copy_row_to(k,Dk);
    A.copy_col_to(k,Ak);
    LD += universal_codelength(M,Dk.weight());
    LA += universal_codelength(N,Ak.weight());
  }
  Dk.destroy();
  Ak.destroy();
  return LE+LD+LA;
}

idx_t learn_model_mdl_forward_selection(binary_matrix& X,
					binary_matrix& E, 
					binary_matrix& D, 
					binary_matrix& A) {
  const idx_t M = E.get_cols();
  const idx_t N = E.get_rows();
  idx_t K = D.get_rows();
  learn_model_inner(X,E,D,A);
  binary_matrix nextAtom(1,M);
  binary_matrix nextCoefs(N,1);
  binary_matrix currD(D),currA(A),currE(E);  
  binary_matrix nextD,nextA,nextE(N,M);
  idx_t bestK = K;
  idx_t bestL = model_codelength(E,D,A);
  idx_t stuck = 0;
  idx_t sumStuck = 0;
  idx_t allStuck = 0;
  do {
    //
    // initialize curr atom and associated coefs.
    //
    idx_t currL = model_codelength(currE,currD,currA);
    int dif = int(currL) - int(bestL);
    int dev = allStuck > 0 ? (sumStuck/allStuck) : 0;
    std::cout << "currK=" << K << " currL=" << currL << " bestK=" << bestK << " bestL=" << bestL << " stuck=" << stuck << " dif=" << dif << " dev=" << dev << std::endl;
    initialize_model(currE,nextAtom,nextCoefs);
    //learn_model_inner(currE,nextE,nextAtom,nextCoefs);
    //std::cout << nextAtom << std::endl;
    //std::cout << nextAtom.weight() << std::endl;
    //
    // create curr dictionary with old atoms plus new one
    // idem with coefs
    //
    // this is painful
    nextD.destroy();
    nextD.allocate(K+1,M);
    nextD.set_submatrix(0,0,currD);
    nextD.set_submatrix(K,0,nextAtom);
    currD.destroy();
    currD.allocate(K+1,M);
    nextD.copy_to(currD);
    nextD.destroy();

    nextA.destroy();
    nextA.allocate(N,K+1);
    nextA.set_submatrix(0,0,currA);
    nextA.set_submatrix(0,K,nextCoefs);
    currA.destroy();
    currA.allocate(N,K+1);
    nextA.copy_to(currA);
    nextA.destroy();

    learn_model_inner(X,currE,currD,currA);
    currL = model_codelength(currE,currD,currA);
    if ((currL + dev) < bestL) {
      stuck = 0;
      bestL = currL;
      D.destroy();
      D.allocate(K+1,M);
      currD.copy_to(D);
      A.destroy();
      A.allocate(N,K+1);
      currA.copy_to(A);
      currE.copy_to(E);
      bestK = K+1;
    } else { 
      stuck++;
      allStuck++;
      sumStuck += (currL-bestL);
      if (stuck >= 10) {
	std::cout << "No further improvement." << std::endl;
	break;
      }
    }
    K++;
    } while (stuck < 10);
  currD.destroy();
  currA.destroy();
  currE.destroy();
  nextE.destroy();
  nextAtom.destroy();
  nextCoefs.destroy();
  return bestL;
}

idx_t learn_model_mdl_backward_selection(binary_matrix& X,
					 binary_matrix& E, 
					 binary_matrix& D, 
					 binary_matrix& A) {
  const idx_t M = E.get_cols();
  const idx_t N = E.get_rows();
  idx_t K = D.get_rows();
  learn_model_inner(X,E,D,A);
  //mul(A,false,D,false,E);
  //add(E,X,E);
  idx_t bestL = model_codelength(E,D,A);
  idx_t bestK = K;
  idx_t currL = bestL;
  binary_matrix Dk(1,M);
  binary_matrix Ak(1,N);
  binary_matrix nextD,currD(D);
  binary_matrix nextA,currA(A);
  binary_matrix nextE(N,M);
  binary_matrix AkDk(N,M);
  
  idx_t stuck = 0;
  idx_t sumStuck = 0;
  idx_t allStuck = 0;
  for (; K > 0; K--) {
    int dif = int(currL) - int(bestL);
    int dev = allStuck > 0 ? (sumStuck/allStuck) : 0;
    std::cout << "currK=" << K << " currL=" << currL << " bestK=" << bestK << " bestL=" << bestL << " stuck=" << stuck << " dif=" << dif << " dev=" << dev << std::endl;
    idx_t nextk = 0;
    idx_t nextL = ~(1UL<<(sizeof(idx_t)-1));
    for (idx_t k = 0; k < K; k++) {
      currD.copy_row_to(k,Dk);
      currA.copy_col_to(k,Ak);
      mul(Ak,true,Dk,false,AkDk);
      // codelength of nex
      add(AkDk,E,nextE); // residual associated to removing Dk
      idx_t tmpL = model_codelength(nextE,currD,currA);
      tmpL -= universal_codelength(M,Dk.weight());
      tmpL -= universal_codelength(N,Ak.weight());
      if (tmpL < nextL) {
	nextL = tmpL;
	nextk = k;
      }
    }
    //
    // create new dictionary and coefficients without discarded atom bestk
    //    
    nextD.destroy();
    nextA.destroy();
    if (K > 1) {
      nextD.allocate(K-1,M);
      nextA.allocate(N,K-1);
      for (idx_t k = 0; k < nextk; k++) {
	currD.copy_row_to(k,Dk);
	nextD.set_row(k,Dk);
	currA.copy_col_to(k,Ak);
	nextA.set_col(k,Ak);
      }
      for (idx_t k = nextk+1; k < K; k++) {
	currD.copy_row_to(k,Dk);
	nextD.set_row(k-1,Dk);
	currA.copy_col_to(k,Ak);
	nextA.set_col(k-1,Ak);
      }
      learn_model_inner(X,nextE,nextD,nextA);
      nextL = model_codelength(nextE,nextD,nextA);
    } else {
      nextL = model_codelength(nextE,nextD,nextA);
    }

    if (nextL + dev < bestL)  {
      // keep best solution so far
      if (K == 1) {
	std::cout << "Resulted in empty model!" << std::endl;
	D.destroy();
	A.destroy();
	X.copy_to(E);
	break;
      }
      stuck = 0;
      bestK = (K-1);
      bestL = nextL;
      D.destroy();
      D.allocate(K-1,M);
      nextD.copy_to(D);
      A.destroy();
      A.allocate(N,K-1);
      nextA.copy_to(A);
      nextE.copy_to(E);
    } else { // not better, do not update solution
      stuck++;
      allStuck++;
      sumStuck += (nextL-bestL);
      if (stuck >= 10) {
	std::cout << "No further improvement." << std::endl;
	break;
      }
    }
    currD.destroy();
    currD.allocate(K-1,M);
    nextD.copy_to(currD);
    
    currA.destroy();
    currA.allocate(N,K-1);
    nextA.copy_to(currA);

    currL = nextL;
  }

  currD.destroy(); nextD.destroy();
  currA.destroy(); nextA.destroy();
  nextE.destroy();
  AkDk.destroy();
  Ak.destroy();
  Dk.destroy();
  return bestL;
}

idx_t learn_model_mdl_full_search(binary_matrix& X,
				  binary_matrix& E, 
				  binary_matrix& D, 
				  binary_matrix& A) {
  const idx_t M = E.get_cols();
  const idx_t N = E.get_rows();
  idx_t K = D.get_rows();
  
  binary_matrix  candE(N,M);
  idx_t bestL = 1UL<<30; // will not work on 16 bit machines
  idx_t bestk = 0;
  for (idx_t k = 20; k <= K; k+=20) {
  //for (idx_t k = 1; k < K; k++) {
    binary_matrix candD(k,M),candA(N,k);
    initialize_model(X,candD,candA);
    learn_model_inner(X,candE,candD,candA);
#define REPS 10
#if 1
    idx_t aux[REPS];
    std::cout << "K=" << k;
    for (idx_t I = 0; I < REPS; I++) {
      random_seed = (random_seed*31 ) % 17;
      initialize_model(X,candD,candA);
      learn_model_inner(X,candE,candD,candA);
      aux[I] = model_codelength(candE,candD,candA);
      std::cout << " " << aux[I];
    }
    idx_t candL = *std::min_element(aux+0,aux+REPS);
    std::cout << " " << candL << std::endl;
#else
    idx_t candL = model_codelength(candE,candD,candA);
    std::cout << "K=" << k << " candL=" << candL << std::endl;
#endif
    if (candL < bestL) {
      bestL = candL;
      bestk = k;
      candE.copy_to(E);

      D.destroy();
      D.allocate(k,M);
      candD.copy_to(D);

      A.destroy();
      A.allocate(N,k);
      candA.copy_to(A);
    }
    candD.destroy();
    candA.destroy();
  }
  candE.destroy();
  std::cout << "bestL=" << bestL << " bestk=" << bestk << std::endl;
  return bestL;
}
