#include "pbm.h"
#include <cstdio>

int main(int argc, char** argv) {
  idx_t rows,cols;
  FILE* fimg;
  fimg = fopen("data/camera.pbm","r");
  if (!fimg) return -1;
  read_pbm_header(fimg,rows,cols);
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;
  binary_matrix A(rows,cols);
  read_pbm_data(fimg,A);
  fclose(fimg);
  fimg = fopen("copy.pbm","w");
  if (!fimg) return -2;
  write_pbm(A,fimg);
  fclose(fimg);
  std::cout << std::endl << A << std::endl;
}
