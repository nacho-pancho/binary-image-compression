#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
//#include <gsl.h>
#include "gsl/gsl_sf_gamma.h"
#include "GolombCoder.h"
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
/*
 * Allows for choosing between predictive and non-predictive coding of residual
 */



int main(int argc, char **argv) {   
  idx_t rows,cols;
  int res;
  FILE* fimg;
  fimg = fopen(argc > 1 ? argv[1]: "data/test.pbm","r");
  idx_t W = 16;
  set_grid_width(W);
  if (!fimg) return -1;
  res = read_pbm_header(fimg,rows,cols);
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;
  binary_matrix I(rows,cols);
  read_pbm_data(fimg,I);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1); }
  fclose(fimg);
  idx_t Ny = (W-1+rows)/W;
  idx_t Nx = (W-1+cols)/W;
  binary_matrix P;
  for (idx_t i = 0; i < Ny; i++) {
    for (idx_t j = 0; j < Nx; j++) {
      const int i0 = i*W;
      const int j0 = j*W;
      P = I.get_submatrix(i0,i0+W,j0,j0+W);
      I.set_submatrix(i0,j0,P); // debug!
    } // main loop: columns
  } // main loop: rows
  fimg = fopen("diff.pbm","w");
  write_pbm(I,fimg);
  fclose(fimg);
  I.destroy();
  return 0;
}
