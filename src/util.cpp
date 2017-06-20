#include "binmat.h"
#include "pbm.h"
#include <cmath>
#include <cstdio>
#include "util.h"

void counting_sort(aux_t* s, idx_t n) {  
  idx_t maxs = 0; 

  //  std::cout << "BEFORE SORT:";
  for (idx_t i = 0; i < n; i++) {
    if (s[i].first > maxs)
      maxs = s[i].first;
    //    std::cout << '(' << s[i].first << ' ' << s[i].second << "), ";
  } 
  //  std::cout << std::endl;
  idx_t count[maxs+2];
  aux_t scratch[n];
  std::fill(count,count+maxs+2,0);
  
  /* Set the value of count[i] to the number of
   * elements in array with value i+min-1. */
  for(idx_t i=0; i<n; i++) {
    count[s[i].first+1]++;
  }

  /* Update count[i] to be the number of
   * elements with value less than i+min. */
  for(idx_t i=1; i <= maxs+1; i++)
    count[i] += count[i-1];

  /* Copy the elements of array into scratch in
   * stable sorted order. */
  for(int i=n-1; i>=0; i--) {
    int c = s[i].first;
    int j = s[i].second;
    scratch[count[c]].first = c;
    scratch[count[c]].second = j;
    /* Increment count so that the next element
     * with the same value as the current element
     * is placed into its own position in scratch. */
    count[c]++;
  }
  //  std::cout << "AFTER SORT:";
  for(idx_t i=0; i<n; i++) {
    s[i] = scratch[i];
    //    std::cout << '(' << s[i].first << ' ' << s[i].second << "), ";
  }

  //  std::cout<<  std::endl;
}

void render_mosaic(const binary_matrix& D, const char* fname) {
  const idx_t m = D.get_cols();
  const idx_t n = D.get_rows();
  const idx_t w = (idx_t) sqrt(double(m));
  const idx_t gn = (idx_t) ceil(sqrt(double(n)));
  const idx_t gm = ceil(double(n)/double(gn));
  const idx_t gw = w + 1;
  const idx_t im = gm*gw;
  const idx_t in = gn*gw;
  binary_matrix I(im,in);
  binary_matrix V(1,m);
  binary_matrix P(w,w);
  idx_t li = 0;
  I.clear();
  for (idx_t i = 0; i < gm; i++) {
    for (idx_t j = 0; j < gn; j++, li++) {
      if (li >= n) break;
      D.copy_row_to(li,V);
      P.set_vectorized(V);
      I.set_submatrix(gw*i,gw*j,P);
    }
    if (li >= n) break;
  }
  FILE* fimg = fopen(fname,"w");
  write_pbm(I,fimg);
  fclose(fimg);
  I.destroy();
  V.destroy();
  P.destroy();
}

int write_pbm(binary_matrix& A, const char* fname) {
  FILE* fimg = fopen(fname,"w");
  if (!fimg) return -2;
  write_pbm(A,fimg);
  fclose(fimg);
  return 0;
}
