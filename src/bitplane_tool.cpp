    #include <cstdio>
#include "pbm.h"
#include "pnm.h"

const char * fname = "data/img/camera.pgm";

int main(int argc, char** argv) {
  int type,rows,cols,maxval;
  FILE* fimg;
  char ofname[128];
  if (argc > 1) { fname = argv[1]; }
  fimg = fopen(fname,"r");
  if (!fimg) {
        std::cerr << "Error opening " << fname << " for reading." << std::endl;
      return -1;
    }      
  read_pnm_header(fimg,type,cols,rows,maxval);
  std::cout << "rows=" << rows << " cols=" << cols << " maxval=" << maxval << std::endl;
  pixel_t* gray_img = new pixel_t[rows*cols];
  binary_matrix A(rows,cols);

  read_pgm_data(fimg,type,rows,cols,maxval,gray_img);
  fclose(fimg);
  for (int b = idx_t(1), bi = 0; b < maxval; b <<= 1, bi++) {
    std::cout << "b=" << b << std::endl;
    for (int i = 0, li = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++,li++) {
	A.set(i,j,gray_img[li] & b);
      }
    }
    snprintf(ofname,128,"plane_%02d.pbm",bi);
    fimg = fopen(ofname,"w");
    if (!fimg) { 
        std::cerr << "Error opening " << ofname << " for writing." << std::endl;
        return -2;
    }
    write_pbm(A,fimg);
    fclose(fimg);
  }
  
  A.destroy();
  delete[] gray_img;
}
