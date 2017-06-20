#ifndef EG_H
#define EG_H
//#include "BinaryFileWriter.h"
//#include "BinaryFileReader.h"

/* DUMMY Golomb Runlength coder, for coding patch weights */
class EG {
public:
  EG() { g = 1; blockSize = 1; lutIndex = 0; };

  //protected: // for debugging
    unsigned g;        // log2(block size)
    unsigned blockSize;// current block size
    int lutIndex; // current index into LUT
    void incBlockSize();
    void decBlockSize();
};

class EGCoder : public EG {
public:
// EGCoder(BinaryFileWriter *f):EG(),file(f) {}
 EGCoder():EG(),bitcount() {}
  void codeRun(int len, bool eol); 
  unsigned long bitcount;
private:
//    BinaryFileWriter * file;
};

class EGDecoder : public EG {
public:
 //EGDecoder(BinaryFileReader *f):EG(),file(f) { }
 EGDecoder():EG() { }
    int decodeRun(int maxlen);
private:
//    BinaryFileReader * file;
};

#endif
