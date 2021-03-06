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

#define XOR(a,b) ( (!(a) && (b)) ||  ((a) && !(b)) )

void med(const binary_matrix& P, binary_matrix& pP) {
  // can be very quickly implemented at block level using binary operators,
  // but I have to sit down and work the math. It's very easy actually...
  for (idx_t j = 1; j < P.get_cols();++j) {
    pP.set(0,j,XOR(P.get(0,j-1),P.get(0,j))); // first row is order 0 pred
  }
  for (idx_t i = 1; i < P.get_rows();++i) {
    pP.set(i,0,XOR(P.get(i-1,0),P.get(i,0))); // first row is order 0 pred      
    for (idx_t j = 1; j < P.get_cols();++j) {
      pP.set(i,j,XOR(XOR(P.get(i-1,j-1),P.get(i,j-1)),XOR(P.get(i-1,j),P.get(i,j))) ); // first row is order 0 pred      
    }
  }
}

int main(int argc, char **argv) {   
  idx_t rows,cols;
  int res;
  FILE* fimg;
  fimg = fopen(argc > 1 ? argv[1]: "data/test.pbm","r");
  const int W = argc > 2 ? atoi(argv[2]) : 16;
  const int R = argc > 4 ? atoi(argv[4]) : 128;
  // a good threshold is one that ensures that a non-matched coding is cheaper
  // than even the best matched coding
  // since the latter involves describing the index of the best patch, which
  // takes about log2 (2*R+1)*R bit
  // and a non-matched coding with weight w takes about w*log2(M) bits
  // we have that any weight below log2 (2*R+1)*R / log2 (M) is "perfect"
  // 
  const idx_t M = W*W;
  const idx_t goodT = ceil(log2(double((2*R+1)*R)))/ceil(log2(double(M)));
  const idx_t T = argc > 3 ? atoi(argv[3]) : goodT;
  std::cout << "T=" << T << std::endl;
  //  const int O = argc > 5 ? atoi(argv[5]) : 0; // always non-predictive
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
  //std::cout << "D:<<" << D << std::endl;
  //std::cout << "iD:<<" << iD << std::endl;
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
  char map[Ny][Nx];
  GolombCoder golomb_match,golomb_nomatch;
  for (idx_t i = 0; i < Ny; i++) {
    for (idx_t j = 0; j < Nx; j++,li++) {
      const int i0 = i*W;
      const int j0 = j*W;
      int i2;
      P = I.get_submatrix(i0,i0+W,j0,j0+W);

      //I.set_submatrix(i0,j0,P); // debug!
      idx_t besti = 0, bestj = 0, bestd = M+1;
      
      bool bestinv = (P.weight() - M) < (P.weight());
      bool perfect_match = (P.weight() <= T) || (P.weight() >= (M-T));
      const int mini = i0>R ? (i0-R) : 0;
      const int minj = (j0>R) ? (j0-R) : 0;
      const int maxj = ((j0+R)>(int(cols)-W)) ? (cols-W): (j0+R);
      const int mini2 = (i0 > W) ? (i0-W) : 0;
      const int maxj2 = (j0 > W) ? (j0-W) : 0;
      const int search_win_size = (i0-mini2)*(maxj2-minj) + (mini2-mini)*(maxj-minj);
      std::cout << std::endl << mini2 << ' ' << i0 << ' '; 
      std::cout << minj << ' ' << maxj2 << ' '; 
      std::cout << std::endl << mini << ' ' << mini2 << ' '; 
      std::cout << minj << ' ' << maxj << ' '; 
      //
      // patches behind current one on same row
      //
      for (i2 = i0; (i2 >= mini2) && !perfect_match; i2--) {
	for (int j2 = maxj2; j2 >= minj; j2--) {
//	  std::cout << i2 << ' ' << j2 << std::endl;
	  P2 = I.get_submatrix(i2,i2+W,j2,j2+W);
	  idx_t d = dist(P,P2);
	  idx_t inv;
	  if ((M-d) < d) {
	    inv = true;
	    d = M-d;
	  }
	  if (d < bestd) {
	    bestinv = inv;
	    bestd = d;
	    besti = i2;
	    bestj = j2;
	    if (bestd <= T)  {
	      perfect_match = true;
	      break;
	    }
	  }	
	}
      }
      //
      // all patches above
      //
      for (i2 = i0-W; (i2 >= mini) && !perfect_match; i2--) {
	for (int j2 = maxj; j2 >= minj; j2--) {
	  P2 = I.get_submatrix(i2,i2+W,j2,j2+W);
	  // could consider worst as 1/2 instead of M if one could use inverted patches
	  idx_t d = dist(P,P2);
	  idx_t inv;
	  if ((M-d) < d) {
	    inv = true;
	    d = M-d;
	  }
	  if (d < bestd) {
	    bestinv = inv;
	    bestd = d;
	    besti = i2;
	    bestj = j2;
	    if (bestd <= T)  {
	      perfect_match = true;
	      break;
	    }
	  }	
	}
      } 
      std::cout << "\n===========================================================================\n" << std::endl;
      std::cout << "i0=" << (i*W) << " j0=" << (j*W);
      std::cout << " besti=" << besti << " bestj=" << bestj << " bestd=" << bestd << " inv=" << bestinv << std::endl;
      histi[besti]++;
      histj[bestj]++;
      histr[ int(sqrt(double(besti*besti+bestj*bestj))) ]++;
      //std::cout << "patch: " << P << std::endl;
      //      std::cout << "best match:" << P2 << std::endl;
      // we found it best to invert the patch
      if (bestinv) {
	P.flip();
      }
      if (bestd <= idx_t(M)) {
	P2 = I.get_submatrix(besti,besti+W,bestj,bestj+W);
	add(P,P2,P3);
      } else {
	P3 = P.get_copy();
      }
      const idx_t match_nonpred_weight = P3.weight();
      const idx_t nomatch_nonpred_weight = P.weight();

      // encode P3 differentially
#if 0 
      binary_matrix V = P3.get_vectorized();
      binary_matrix dP3(1,M);
      mul(V,false,D,false,dP3);
      binary_matrix dP(1,M);
      V = P.get_vectorized();
      mul(V,false,D,false,dP);
#else // MED
      binary_matrix dP(W,W);
      binary_matrix dP3(W,W);
      med(P,dP);
      med(P3,dP3);
#endif
      const idx_t match_pred_weight = dP3.weight();
      const idx_t nomatch_pred_weight = dP.weight();

      const idx_t idx_len = ceil(log2(search_win_size));
      std::cout << "Search cost: size=" << search_win_size << " len=" << idx_len << std::endl;
      //      std::cout << "w(P3)=" << P3.weight() << " w(dP3)=" << dP3.weight() << std::endl;
      std::cout << "weight: nonmatch/nonpred=" << nomatch_nonpred_weight << '\t' ;
      std::cout << "nonmatch/pred=" << nomatch_pred_weight << '\t' ;
      std::cout << "match/nonpred=" << match_nonpred_weight << '\t' ;
      std::cout << "match/pred=" << match_pred_weight << std::endl;
      
      //std::cout << "difference:" << P3 << std::endl;
      // 2 fixed bits: one for indicating match/nomatch, other for pred/nonpred, 
      const idx_t nomatch_nonpred_len = 2 + enumL(M,nomatch_nonpred_weight);
      const idx_t nomatch_pred_len = 2 + enumL(M,nomatch_pred_weight);
      // in the case of match, there is another for invert/noinvert
      const idx_t match_nonpred_len = 3 + idx_len + enumL(M,match_nonpred_weight);
      const idx_t match_pred_len = 3 + idx_len + enumL(M,match_pred_weight);
      std::cout << "len: nonmatch/nonpred=" << nomatch_nonpred_len << '\t' ;
      std::cout << "nonmatch/pred=" << nomatch_pred_len << '\t' ;
      std::cout << "match/nonpred=" << match_nonpred_len << '\t' ;
      std::cout << "match/pred=" << match_pred_len << std::endl;

#define MIN(a,b) ((a) < (b) ? (a) : (b))

      idx_t match_len, nomatch_len;
      idx_t match_weight, nomatch_weight;
      binary_matrix residual_match(W,W), residual_nomatch(W,W);
      char match_mode, nomatch_mode;

      if (match_nonpred_len > match_pred_len) {
	match_len = match_pred_len;
	match_weight = match_pred_weight;
	residual_match = dP3;
	match_mode = 'X';
      } else {
	match_len = match_nonpred_len;
	match_weight = match_nonpred_weight;
	residual_match = P3;
	match_mode = 'x';
      }

      if (nomatch_nonpred_len > nomatch_pred_len) {
	nomatch_len = nomatch_pred_len;
	nomatch_weight = nomatch_pred_weight;
	residual_nomatch = dP;
	nomatch_mode = 'O';
      } else {
	nomatch_len = nomatch_nonpred_len;
	nomatch_weight = nomatch_nonpred_weight;
	residual_nomatch = P;
	nomatch_mode = 'o';
      }

      if (nomatch_len > match_len) {
	golomb_match.codeSample(match_weight);
	hist[bestd]++;
	average_weight += match_weight; 
	matches++;
 	L += match_len;
        I.set_submatrix(i0,j0,residual_match);
	std::cout << "mode="<< match_mode << std::endl;
	map[i][j] = match_mode;
      } else {
	golomb_nomatch.codeSample(nomatch_weight);
	L += nomatch_len;
        I.set_submatrix(i0,j0,residual_nomatch);
	std::cout << "mode="<< nomatch_mode << std::endl;
	map[i][j] = nomatch_mode;
      }
    } // main loop: columns
  } // main loop: rows
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
  print_hist(hist,M/4,true);
  //  print_hist(histi,rows,false);
  //print_hist(histj,cols,false);;
  //print_hist(histr,ninf,false);
  std::cout << "MAP:" << std::endl;
  for (idx_t i = 0; i < Ny; i++) {
    for (idx_t j = 0; j < Nx; j++) {
      std::cout << map[i][j];
    }
    std::cout<<std::endl;
  }

  fimg = fopen("diff.pbm","w");
  write_pbm(I,fimg);


  P.destroy();
  P2.destroy();
  P3.destroy();
  I.destroy();
  I2.destroy();
  return 0;
}
