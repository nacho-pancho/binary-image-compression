#include <cstdio>
#include "pbm.h"
#include "pnm.h"
#include <iomanip>

const char * ifbasename = "data/camera_%02d.pbm"; 
const char * ofname = "data/camera_rec.pgm"; 

int main(int argc, char** argv) {
  idx_t type,rows,cols;
  FILE* fimg;
  char ifname[128];
  if (argc > 1) { ifbasename = argv[1]; }
  if (argc > 2) { ofname = argv[2]; }

  unsigned plane = 0;
  snprintf(ifname,128,ifbasename,plane);
  fimg = fopen(ifname,"r");
  if (!fimg) {
      std::cerr << "Error reading  image " << ifname <<std::endl;
      return -1;
  }
  type = 5;
  read_pbm_header(fimg,rows,cols);
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;
  pixel_t* gray_img = new pixel_t[rows*cols];
  std::fill(gray_img,gray_img+rows*cols,0);
  binary_matrix A(rows,cols);
  A.clear();
  read_pbm_data(fimg,A);
  fclose(fimg);
  unsigned mask = 0x01;
  do {
    std::cout << "plane=" << std::dec << plane << " mask=" << std::hex << mask << std::endl;
    for (unsigned i = 0, li = 0; i < rows; i++) {
      for (unsigned j = 0; j < cols; j++,li++) {
          if (i == 0 && j < 32) std::cout << A.get(i,j) << ',';
        if (A.get(i,j))
            gray_img[li] |= mask;
      }
    }
    plane++;
    mask <<= 1;
    snprintf(ifname,128,ifbasename,plane);  
    fimg = fopen(ifname,"r");
    if (fimg) {
        read_pbm_data(fimg,A);
        fclose(fimg);
    }
  } while (fimg);
  A.destroy();
  write_pgm(gray_img,type,cols,rows,mask,ofname);
  delete[] gray_img;
}
