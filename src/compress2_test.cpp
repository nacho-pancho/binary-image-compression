#include "binmat.h"
#include "pbm.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

//#include <gsl.h>
#include "gsl/gsl_sf_gamma.h"
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
/*
COMPRESSION ALGORITHMS:
non-dictionary: send index of best patch + encoded difference using enumerative coding
if difference is too large, use no other patch
+ no parameters
- index of best patch can be quite large
- very slow
dictionary: 
if diff. with best patch in dictionary above threshold, add to dictionary and encode original
send index in dictionary plus difference encoded using enumerative
parameter: threshold

V2: dictionary based
 */

#define COSMOS_2E	5.436563656918090181591196596855297684669
#define COSMOS_2PI	6.283185307179586231995926937088370323181
#define COSMOS_2EPI	17.07946844534713193297648103907704353333
#define COSMOS_LOG2E	1.442695040888963387004650940070860087872
#define COSMOS_LOG2PI	1.651496129472318719066947778628673404455
#define COSMOS_LOG2EPI	4.0941911703612818840269937936682254076

//#define W 7 // diabolic

double enumL(const idx_t n,
	     const idx_t r) { // number of nonzeroes
  return r>0 ? gsl_sf_lnchoose ( n,r ) *COSMOS_LOG2E : 0.0 ;
}

typedef struct { idx_t i,j,hits; } coord_t;

int main(int argc, char **argv) {   
  idx_t rows,cols;
  int res;
  FILE* fimg;
  fimg = fopen(argc > 1 ? argv[1]: "data/einstein.pbm","r");
  idx_t W = argc > 2 ? atoi(argv[2]) : 16;

  idx_t hist[W*W];
  for (idx_t i = 0; i < (W*W) ; i++)
    hist[i] = 0;

  std::vector<coord_t> dict;

  set_grid_width(W);
  if (!fimg) return -1;
  res = read_pbm_header(fimg,rows,cols);
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;
  binary_matrix I(rows,cols);

  read_pbm_data(fimg,I);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1); }
  fclose(fimg);
  //std::cout << "==== ORIGINAL =====\n" << std::endl;
  //std::cout << std::endl << I << std::endl;
  binary_matrix P,P2,I2;
  idx_t Ny = (W-1+rows)/W;
  idx_t Nx = (W-1+cols)/W;
  //  std::cout << "Nx=" << Nx << " Ny=" << Ny << std::endl;
  //  binary_matrix D(Nx*Ny,M);
  idx_t li = 0;
  I2 = I.get_copy();
  binary_matrix P3(W,W);
  double L = 0;
  for (idx_t i = 0; i < Ny; i++) {
    for (idx_t j = 0; j < Nx; j++,li++) {
      const idx_t i0 = i*W;
      const idx_t j0 = j*W;
      P = I.get_submatrix(i0,i0+W,j0,j0+W);
      idx_t bestk = 0, bestd = W*W;
      bool perfect_match = false;
      for (idx_t k = 0; (k < dict.size()) && !perfect_match; k++) {
	coord_t c = dict[k];
	P2 = I.get_submatrix(c.i,c.i+W,c.j,c.j+W);
	const idx_t d = dist(P,P2);
	if (d < bestd) {
	  bestd = d;
	  bestk = k;
	}	
	if (bestd == 0)  {
	  perfect_match = true;
	  break;
	}
      }
      std::cout << "\n===========================================================================\n" << std::endl;
      std::cout << "i=" << (i*W) << " j=" << (j*W);
      std::cout << " bestk=" << bestk << " bestd=" << bestd << " |D|=" << dict.size() << " log2|D|=" << ceil(log2(double(dict.size()))) << std::endl;
//      std::cout << "patch: " << P << std::endl;
      const idx_t nomatch_len = 1 + enumL(W*W,P.weight()) +  0.5*log2 ( double(W*W) );
      if (dict.size()== 0) {
	const coord_t newc = {i,j,0};
	dict.push_back(newc);
	L += nomatch_len;
	continue;
      } 
      coord_t c = dict[bestk];
      P2 = I.get_submatrix(c.i,c.i+W,c.j,c.j+W);
 //     std::cout << "best match:" << P2 << std::endl;
      add(P,P2,P3);
 //     std::cout << "difference:" << P3 << std::endl;
      const idx_t match_len = 1 + ceil(log2(dict.size())) + enumL(W*W,bestd) + 0.5*log2 (double(W*W));
      // we can *vastly* improve on the last term 0.5log2 if histogram of weight is used
      std::cout << "nomatch len=" << nomatch_len << " match_len=" << match_len;
      if (nomatch_len > match_len) {
	std::cout << " USE MATCH!" << std::endl;
	hist[bestd]++;
	L += match_len;
      } else { // no match
	std::cout << " ADD TO DICT!" << std::endl;
	const coord_t newc = {i,j,0};
	dict.push_back(newc);
	L += nomatch_len;
      }
    }    
  }
  std::cout << " COMP CODELENGTH: " << L << std::endl;
  std::cout << " RAW CODELENGTH: " << rows*cols << std::endl;
  std::cout << "RATIO: " << (100.0*L/(rows*cols)) << std::endl;
  for (idx_t i = 0; i < (W*W) ; i++) {
    std::cout << i << ":";
    for (idx_t j = 0; j < hist[i] ; j++) {
      std::cout << "#";
    }
    std::cout << std::endl;
  }
  P.destroy();
  P2.destroy();
  I.destroy();
  I2.destroy();
  return 0;
}
