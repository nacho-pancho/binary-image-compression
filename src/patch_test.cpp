#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>

//#define W 7 // diabolic

int main(int argc, char **argv) {   
  idx_t rows,cols;
  int res;
  std::cout << "ONES=" << bm_bitset(ONES) << "ZEROES=" << bm_bitset(ZEROES) << std::endl;
  std::cout << "MSB=" << bm_bitset(MSB) << " LSB=" << bm_bitset(LSB) << std::endl;
  FILE* fimg;
  fimg = fopen(argc > 1 ? argv[1]: "data/test.pbm","r");
  idx_t W = argc > 2 ? atoi(argv[2]) : 16;
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
  std::cout << "Nx=" << Nx << " Ny=" << Ny << std::endl;
  binary_matrix D(Nx*Ny,M);
  idx_t li = 0;
  binary_matrix I2(rows,cols);
  I2.clear();
  for (idx_t i = 0; i < Ny; i++) {
    for (idx_t j = 0; j < Nx; j++,li++) {
      //     std::cout << "n=" << li << std::endl;
      P = I.get_submatrix(i*W,(i+1)*W,j*W,(j+1)*W);
      P2 = P.get_copy();
            std::cout << P << std::endl;
      V = P.get_vectorized();
            std::cout <<  V<< std::endl;
      P.set_vectorized(V);
            std::cout << P << std::endl;
      if ((i == 1) && (j >= 6)) { 
	// for debugging
      }
      if (dist(P,P2)> 0) {
	std::cout << "Difference after vect! " << dist(P,P2) << std::endl;
      }
      I2.set_submatrix(i*W,j*W,P);
      D.set_row(li,V);
    }
  }
  std::cout << "==== REWRITTEN =====\n" << std::endl;
  std::cout << std::endl << I2 << std::endl;
  fimg = fopen("rewritten.pbm","w");
  if (!fimg) return -2;
  write_pbm(I,fimg);
  fclose(fimg);
  if (dist(I,I2)) {
    std::cout << "DIFFER after first pass!" << std::endl;
  }
  std::cout << "DICTIONARY:\n" << D << std::endl;
  I.clear();
  li = 0;
  for (idx_t i = 0; i < Ny; i++) {
    for (idx_t j = 0; j < Nx; j++,li++) {
      //      std::cout << "n=" << li << std::endl;
      V = D.get_row(li);
      P = binary_matrix(W,W);
      P.set_vectorized(V);
      //      std::cout << P << std::endl;
      I.set_submatrix(i*W,j*W,P);
    }
  }
 
  std::cout << "==== RECONSTRUCTED =====\n" << std::endl;
  std::cout << std::endl << I << std::endl;
  fimg = fopen("reconstructed.pbm","w");
  if (!fimg) return -2;
  write_pbm(I,fimg);
  fclose(fimg);
  if (dist(I,I2)) {
    std::cout << "DIFFER after first pass!" << std::endl;
  }
  V.destroy();
  D.destroy();
  P.destroy();
  P2.destroy();
  I.destroy();
  I2.destroy();
  return 0;
}
