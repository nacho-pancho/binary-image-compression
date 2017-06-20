#include "eg.h"
static const short EGLUT[32] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 9, 10, 11, 12, 13, 14, 15};

void EG::incBlockSize() {
  if (lutIndex < 32)
    lutIndex++;
  g = EGLUT[lutIndex];
  blockSize = 1UL << g;
  //printf("inc: idx=%d g=%d bs=%d\n",lutIndex,g,blockSize);
}

void EG::decBlockSize() {
  if (lutIndex > 0)
    lutIndex--;
  g = EGLUT[lutIndex];
  blockSize = 1UL << g;
  //printf("dec: idx=%d g=%d bs=%d\n",lutIndex,g,blockSize);
}

void EGCoder::codeRun(int len, bool eol) { 
  //printf("coding run of length %d (eol=%d)",len,eol);
  while (len >= blockSize) {
    len -= blockSize;
//    file->writeBits(1,1); 
//    incBlockSize();
    bitcount++;
  }
  if (eol) { // run was broken by EOL
    //file->writeBits(1,1);
    bitcount++;
  } else { // len < blocksize
//    file->writeBits(0,1);
//    file->writeBits(len,g);
    bitcount+=(g+1);
    decBlockSize();
  }
}

#if 0

int EGDecoder::decodeRun(int maxlen) {
  int len = 0;
  int fullblock;
  while ( fullblock = file->readBits(1) ) { 
    len += blockSize;
    if (len > maxlen) {
      return maxlen;
    }
    incBlockSize();
  }
  //printf("Decoded run of length %d (eol=%d)\n",len + file->readBits(g),false);
  len += file->readBits(g);
  decBlockSize();
  return len;
}
#endif
