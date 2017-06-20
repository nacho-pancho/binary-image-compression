//
//  GolombCoder.cpp
//  compresoreeg
//
//  Created by Ignacio Capurro on 12/13/13.
//  Copyright (c) 2013 Ignacio Capurro. All rights reserved.
//

#include "GolombCoder.h"
#include <cassert>
#include <iostream>

void GolombCoder::binaryEncode(unsigned sample, unsigned k){
  assert(k < sizeof(int)*8);
  if (k > 100) {
    std::cout << "k " << k << endl;
  }
  //  unsigned binary = sample & ~(((unsigned)(-1)) << k); // -1 es la máscara 0xFFFFFFFF
  unsigned unary = sample >> k;
  
  // write the binary part
  //file->writeBits(binary, k);
  // write unary part
  //file->writeZeros(unary);
  //file->writeBits(1, 1);
  bitcount += k + unary + 1;
}

void GolombCoder::codeSample(unsigned sample) {
  binaryEncode(sample,k);    
  samples++;
  accumulatedError += sample;    
  for(k=0 ; (samples<<k) < accumulatedError; k++ );
}
