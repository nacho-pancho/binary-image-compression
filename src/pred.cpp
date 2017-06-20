#include "binmat.h"

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
