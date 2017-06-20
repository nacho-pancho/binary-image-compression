#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>

//#define W 7 // diabolic

int main(int argc, char **argv) {   
  idx_t rows,cols;
  int res;
  FILE* fimg;
  fimg = fopen(argc > 1 ? argv[1]: "data/test.pbm","r");
  idx_t W = argc > 2 ? atoi(argv[2]) : 5;
  set_grid_width(W);
  if (!fimg) return -1;
  res = read_pbm_header(fimg,rows,cols);
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;
  binary_matrix I(rows,cols);
  read_pbm_data(fimg,I);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1); }
  fclose(fimg);
  std::cout << "==== ORIGINAL =====\n" << std::endl;
  std::cout << std::endl << I << std::endl;
  binary_matrix P,V,P2;
  idx_t Ny = (W-1+rows)/W;
  idx_t Nx = (W-1+cols)/W;
  idx_t M = W*W;
  //  std::cout << "Nx=" << Nx << " Ny=" << Ny << std::endl;
  //  binary_matrix D(Nx*Ny,M);
  idx_t li = 0;
  I2 = I.copy();
  binary_matrix P2;
  for (idx_t i = 0; i < Ny; i++) {
    for (idx_t j = 0; j < Nx; j++,li++) {
      const idx_t i0 = i*W;
      const idx_t j0 = j*W;
      P = I.get_submatrix(i0,i1+W,j0,j1+W);
      std::cout << "i=" << i << "j=" << j << std::endl;
      std::cout << P << std::endl;
      idx_t li2 = 0;
      idx_t besti = 0, bestj = 0, bestd = W*W;
      idx_t i2;
      for (i2 = 0; i2 < (i0-W); i++) {
	for (idx_t j2 = 0; j2 < cols; j++,li2++) {
	  P2 = I.get_submatrix(i2,i2+W,j2,j2+W);
	  // could consider worst as 1/2 instead of W*W if one could use inverted patches
	  if (dist(P,P2) < bestd) {
	    bestd = d;
	    besti = i2;
	    bestj = j2;
	  }
	}
      }   
      for (idx_t j2 = 0; j2 < (j0-W); j++,li2++) {
	P2 = I.get_submatrix(i2,i2+W,j2,j2+W);
	if (dist(P,P2) < bestd) {
	  bestd = d;
	  besti = i2;
	  bestj = j2;
	}	
      }
    }
  }
  V.destroy();
  D.destroy();
  P.destroy();
  I.destroy();
  I2.destroy();
  return 0;
}
