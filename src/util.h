#ifndef UTIL_H
#define UTIL_H
#include "binmat.h"

typedef std::pair<idx_t,idx_t> aux_t;
void render_mosaic(const  binary_matrix& D, const char* fname);
void counting_sort(aux_t* s, idx_t n);
int write_pbm(binary_matrix& A, const char* fname);

#endif

