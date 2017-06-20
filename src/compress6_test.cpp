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

- Differs from compress5 in that rsiduals are encoded predictively
 */

#define COSMOS_2E	5.436563656918090181591196596855297684669
#define COSMOS_2PI	6.283185307179586231995926937088370323181
#define COSMOS_2EPI	17.07946844534713193297648103907704353333
#define COSMOS_LOG2E	1.442695040888963387004650940070860087872
#define COSMOS_LOG2PI	1.651496129472318719066947778628673404455
#define COSMOS_LOG2EPI	4.0941911703612818840269937936682254076

double enumL(const idx_t n,
	     const idx_t r) { // number of nonzeroes
  return r>0 ? gsl_sf_lnchoose ( n,r ) *COSMOS_LOG2E : 0.0;
}


void print_hist(idx_t hist[], idx_t n, bool logscale) {
  for (idx_t i = 0; i < n ; i++) {
    std::cout << i << ":";
    const idx_t top = logscale ? ceil(log2(hist[i]+1.)) : hist[i];
    for (idx_t j = 0; j < top ; j++) {
      std::cout << "#";
    }
    std::cout << std::endl;
  }
}

int main(int argc, char **argv) {   
  idx_t rows,cols;
  int res;
  FILE* fimg;
  fimg = fopen(argc > 1 ? argv[1]: "data/test.pbm","r");
  const int W = argc > 2 ? atoi(argv[2]) : 16;
  const idx_t T = argc > 3 ? atoi(argv[3]) : 0;
  const int R = argc > 4 ? atoi(argv[4]) : 10000;
  const idx_t M = W*W;
  idx_t hist[M];
  std::fill(hist,hist+M,0);
  // predictor and inverse predictor matrices
  binary_matrix D(M,M), iD(M,M); 
  D.clear();
  iD.clear();
  for (idx_t i = 0; i < M; i++) {
    D.set(i,i);
    if (i>0) {
      D.set(i-1,i);
    }
    for (idx_t j = i; j < M; j++) {
      iD.set(i,j);
    }
  }
  std::cout << "D:<<" << D << std::endl;
  std::cout << "iD:<<" << iD << std::endl;
  // check invertibility
  //binary_matrix A(M,M);
  //mul(D,false,iD,false,A);
  //std::cout << "D*iD:<<" << A << std::endl;

  set_grid_width(W);
  if (!fimg) return -1;
  res = read_pbm_header(fimg,rows,cols);
  std::cout << "rows=" << rows << " cols=" << cols << std::endl;
  binary_matrix I(rows,cols);
  read_pbm_data(fimg,I);
  if (res !=PBM_OK) { std::cerr << "Error " << res << " reading image."  << std::endl; std::exit(1); }
  fclose(fimg);
  //  std::cout << "==== ORIGINAL =====\n" << std::endl;
  //  std::cout << std::endl << I << std::endl;
  binary_matrix P,P2,I2;
  idx_t Ny = (W-1+rows)/W;
  idx_t Nx = (W-1+cols)/W;
  idx_t histi[rows];
  idx_t histj[cols];
  idx_t ninf = cols < rows ? rows: cols;
  idx_t histr[ninf];
  std::fill(histi,histi+rows,0);
  std::fill(histj,histj+cols,0);
  std::fill(histr,histr+ninf,0);
  idx_t li = 0;
  I2 = I.get_copy();
  binary_matrix P3(W,W);
  double L = 0;
  idx_t average_weight = 0;
  idx_t matches = 0;
  GolombCoder golomb_match,golomb_nomatch;
  for (idx_t i = 0; i < Ny; i++) {
    for (idx_t j = 0; j < Nx; j++,li++) {
      const int i0 = i*W;
      const int j0 = j*W;
      P = I.get_submatrix(i0,i0+W,j0,j0+W);
      idx_t besti = 0, bestj = 0, bestd = M+1;
      int i2;
      bool perfect_match = false;
      const int mini = i0>R ? (i0-R) : 0;
      const int mini2 = (i0 > W) ? (i0-W) : 0;
      const int minj = (j0>R) ? (j0-R) : 0;
      const int maxj = ((j0+R)>(int(cols)-W)) ? (cols-W): (j0+R);
      std::cout << std::endl << mini2 << ' ' << i0 << ' '; 
      std::cout << minj << ' ' << (j0-W) << ' '; 
      for (i2 = i0; (i2 >= mini2) && !perfect_match; i2--) {
	for (int j2 = int(j0-W); j2 >= minj; j2--) {
//	  std::cout << i2 << ' ' << j2 << std::endl;
	  P2 = I.get_submatrix(i2,i2+W,j2,j2+W);
	  const idx_t d = dist(P,P2);
	  if (d < bestd) {
	  //if (d < bestd) {
	    bestd = d;
	    besti = i2;
	    bestj = j2;
	  }	
	  if (bestd <= T)  {
	    perfect_match = true;
	    break;
	  }
	}
      }
      for (i2 = i0-W; (i2 >= mini) && !perfect_match; i2--) {
	for (int j2 = maxj; j2 >= minj; j2--) {
	  P2 = I.get_submatrix(i2,i2+W,j2,j2+W);
	  // could consider worst as 1/2 instead of M if one could use inverted patches
	  const idx_t d = dist(P,P2);
	  if (d < bestd) {
	    bestd = d;
	    besti = i2;
	    bestj = j2;
	  }
	  if (bestd <= T)  {
	    perfect_match = true;
	    break;
	  }
	}
      } 
      std::cout << "\n===========================================================================\n" << std::endl;
      std::cout << "i=" << (i*W) << " j=" << (j*W);
      std::cout << " besti=" << besti << " bestj=" << bestj << " bestd=" << bestd << std::endl;
      histi[besti]++;
      histj[bestj]++;
      histr[ int(sqrt(double(besti*besti+bestj*bestj))) ]++;
      //std::cout << "patch: " << P << std::endl;
      //      std::cout << "best match:" << P2 << std::endl;
      if (bestd <= idx_t(M)) {
	P2 = I.get_submatrix(besti,besti+W,bestj,bestj+W);
	add(P,P2,P3);
      } else {
	P3 = P.get_copy();
      }
#if 0
      // encode P3 differentially
      binary_matrix V = P3.get_vectorized();
      binary_matrix dP3(1,M);
      mul(V,false,D,false,dP3);
      const idx_t match_weight = dP3.weight();
      binary_matrix dP(1,M);
      V = P.get_vectorized();
      mul(V,false,D,false,dP);
      const idx_t nomatch_weight = dP.weight();
#else
      const idx_t match_weight = P3.weight();
      const idx_t nomatch_weight = P.weight();
#endif

      const idx_t idx_len = ceil(log2(li));
      //      std::cout << "w(P3)=" << P3.weight() << " w(dP3)=" << dP3.weight() << std::endl;
      std::cout << "nmweight=" << nomatch_weight << ' ' ;
      std::cout << "mweight=" << match_weight << std::endl;
      
      //std::cout << "difference:" << P3 << std::endl;
      const idx_t nomatch_len = 1 + enumL(M,nomatch_weight);
      const idx_t match_len = 1 + idx_len + enumL(M,match_weight);
      std::cout << "nomatch len=" << nomatch_len << " match_len=" << match_len;
      if (nomatch_len > match_len) {
	golomb_match.codeSample(match_weight);
	std::cout << " USE MATCH!" << std::endl;
	hist[bestd]++;
	average_weight += match_weight; 
	matches++;
 	L += match_len;
        I.set_submatrix(i0,j0,P3);
      } else {
	golomb_nomatch.codeSample(nomatch_weight);
	L += nomatch_len;
      }
    }    
  }
  if (matches == 0) matches++;
  average_weight /= matches;
  std::cout << "\nMATCHES: " << matches  << std::endl;
  std::cout << "AVG. WEIGHT: " << average_weight << std::endl;
  std::cout << "Avg. Golomb/Match: " << golomb_match.bitcount/matches  << std::endl;
  std::cout << "Avg. Golomb/NoMatch: " << golomb_nomatch.bitcount/(Nx*Ny-matches)  << std::endl;
  std::cout << "AVG. WEIGHT: " << average_weight << std::endl;
  L = L+golomb_match.bitcount + golomb_nomatch.bitcount;
  std::cout << "COMP CODELENGTH (bytes): " <<L/8 << std::endl;
  std::cout << "RAW CODELENGTH (bytes): " << (rows*cols)/8 << std::endl;
  std::cout << "RATIO: " << (100.0*L/(rows*cols)) << std::endl;
  print_hist(hist,M,true);
  print_hist(histi,rows,false);
  print_hist(histj,cols,false);;
  print_hist(histr,ninf,false);
  fimg = fopen("diff.pbm","w");
  write_pbm(I,fimg);
  P.destroy();
  P2.destroy();
  P3.destroy();
  I.destroy();
  I2.destroy();
  return 0;
}
